#pragma once

#include <DCA/Utils/FD.h>
#include <DCA/Autodiff/AD_PlaneToSphere.h>
#include <DCA/Interactions/SphereVsSphere.h>

namespace DCA {

class PlaneVsSphere {
public:
    /**
     * @param P Plane pos, Plane normal, Sphere Pos
     * @param props Sphere radius
     */
    static double compute_D(const Vector9d& P, const Vector1d& props) {
        Vector2d radius_padded;
        radius_padded << 0, props;

        Vector3d spherePos = P.segment(6, 3);
        Vector3d planePos = getProjectionOfPoint(P.segment(0, 6), spherePos);
        Vector6d bothPos;
        bothPos << planePos, spherePos;
        return SphereVsSphere::compute_D(bothPos, radius_padded);
    }

    FD_CHECK_dDdP(9, 9, 1, 0, "PlaneVsSphere - dDdP_9");
    // full
    static void compute_dDdP(Vector9d& dDdP, const Vector9d& P,
                             const Vector1d& props) {
        PlaneToSphereDistance_CodeGen::AD_PlaneToSphereDistanceGradient(
            P, props, dDdP);
    }

    FD_CHECK_dDdP(6, 9, 1, 0, "PlaneVsSphere - dDdP_6");
    // to plane
    static void compute_dDdP(Vector6d& dDdP, const Vector9d& P,
                             const Vector1d& props) {
        Vector9d dDdP_full;
        PlaneToSphereDistance_CodeGen::AD_PlaneToSphereDistanceGradient(
            P, props, dDdP_full);
        dDdP = dDdP_full.head(6);
    }

    FD_CHECK_dDdP(3, 9, 1, 6, "PlaneVsSphere - dDdP_3");
    // to sphere
    static void compute_dDdP(Vector3d& dDdP, const Vector9d& P,
                             const Vector1d& props) {
        Vector9d dDdP_full;
        PlaneToSphereDistance_CodeGen::AD_PlaneToSphereDistanceGradient(
            P, props, dDdP_full);
        dDdP = dDdP_full.tail(3);
    }

    FD_CHECK_d2DdP2(9, 9, 1, 0, "PlaneVsSphere - d2DdP2_9");
    static void compute_d2DdP2(Matrix9d& d2DdP2, const Vector9d& P,
                               const Vector1d& props) {
        PlaneToSphereDistance_CodeGen::AD_PlaneToSphereDistanceHessian(P, props,
                                                                       d2DdP2);
    }
    FD_CHECK_d2DdP2(6, 9, 1, 0, "PlaneVsSphere - d2DdP2_6");
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Vector9d& P,
                               const Vector1d& props) {
        Matrix9d d2DdP2_full;
        PlaneToSphereDistance_CodeGen::AD_PlaneToSphereDistanceHessian(
            P, props, d2DdP2_full);
        d2DdP2 = d2DdP2_full.block(0, 0, 6, 6);
    }
    FD_CHECK_d2DdP2(3, 9, 1, 6, "PlaneVsSphere - d2DdP2_3");
    static void compute_d2DdP2(Matrix3d& d2DdP2, const Vector9d& P,
                               const Vector1d& props) {
        Matrix9d d2DdP2_full;
        PlaneToSphereDistance_CodeGen::AD_PlaneToSphereDistanceHessian(
            P, props, d2DdP2_full);
        d2DdP2 = d2DdP2_full.block(6, 6, 3, 3);
    }

public:
    /**
     * @param P Plane pos, Plane normal, Sphere Pos
     * @param props Sphere radius
     */
    static Matrix3d getPointProjectionDerivative(const Vector9d& P,
                                                 const Vector1d& props) {
        Matrix3d mat = Matrix3d::Identity();
        const Vector3d normal = P.segment(3, 3);
        mat.row(0) += -normal.transpose() * normal(0);
        mat.row(1) += -normal.transpose() * normal(1);
        mat.row(2) += -normal.transpose() * normal(2);
        return mat;
    }

    /**
     * @param P Plane pos, Plane normal
     * @todo put in PlaneHelper
     */
    static Vector3d getProjectionOfPoint(const Vector6d& P,
                                         const Vector3d& point) {
        return point + -P.segment(3, 3) * getSignedDistanceToPoint(P, point);
    }

    /**
     * @param P Plane pos, Plane normal
     */
    static double getSignedDistanceToPoint(const Vector6d& P,
                                           const Vector3d& point) {
        return (point - P.head(3)).dot(P.segment(3, 3));
    }

    static void getCartesianEquationCoeffs(const Vector6d& P, double& a,
                                           double& b, double& c, double& d) {
        // normal from 3 to 5
        a = P(3);
        b = P(4);
        c = P(5);
        d = -(a * P(0) + b * P(1) + c * P(2));
    }
};

}  // namespace DCA
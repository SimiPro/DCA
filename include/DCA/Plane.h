#ifndef __DCA_PLANE_H__
#define __DCA_PLANE_H__

#include "AD_PlaneToCapsule.h"
#include "AD_PlaneToPlane.h"
#include "AD_PlaneToSphere.h"
#include "FD.h"
#include "SphereVsSphere.h"
#include "utils.h"

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

class PlaneVsCapsule {
public:
    static double compute_D(const Vector12d& P, const Vector1d& props) {
        Vector3d Pt_line, Pt_plane;
        computeClosestPoints(P, Pt_plane, Pt_line);
        return (Pt_line - Pt_plane).norm() - props(0);  // minus caps. radius
    }

    FD_CHECK_dDdP(6, 12, 1, 0, "PlaneVsCapsule - dDdP_6");
    static void compute_dDdP(Vector6d& dDdP, const Vector12d& P,
                             const Vector1d& props) {
        Vector12d dDdP_full;
        PlaneToCapsuleDistance_CodeGen::AD_PlaneToCapsuleDistanceGradient(
            P, props, curveScale(), sigScale(), dDdP_full);
        dDdP = dDdP_full.head(6);
    }

    FD_CHECK_dDdP(12, 12, 1, 0, "PlaneVsCapsule - dDdP_12");
    static void compute_dDdP(Vector12d& dDdP, const Vector12d& P,
                             const Vector1d& props) {
        PlaneToCapsuleDistance_CodeGen::AD_PlaneToCapsuleDistanceGradient(
            P, props, curveScale(), sigScale(), dDdP);
    }

    FD_CHECK_d2DdP2(12, 12, 1, 0, "PlaneVsCapsule - d2DdP2_12");
    static void compute_d2DdP2(Matrix12d& d2DdP2, const Vector12d& P,
                               const Vector1d& props) {
        PlaneToCapsuleDistance_CodeGen::AD_PlaneToCapsuleDistanceHessian(
            P, props, curveScale(), sigScale(), d2DdP2);
    }

    FD_CHECK_d2DdP2(6, 12, 1, 0, "PlaneVsCapsule - d2DdP2_6");
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Vector12d& P,
                               const Vector1d& props) {
        Matrix12d d2DdP2_full;
        PlaneToCapsuleDistance_CodeGen::AD_PlaneToCapsuleDistanceHessian(
            P, props, curveScale(), sigScale(), d2DdP2_full);
        d2DdP2 = d2DdP2_full.block(0, 0, 6, 6);
    }

public:
    /**
     * @param P plane pos, plane normal, caps start, caps end
     */
    static void computeClosestPoints(const Vector12d& P, Vector3d& Pt_plane,
                                     Vector3d& Pt_line) {
        Vector3d P1 = P.segment(6, 3);  // caps start
        Vector3d P2 = P.segment(9, 3);  // caps end
        Vector3d P1_proj =
            PlaneVsSphere::getProjectionOfPoint(P.segment(0, 6), P1);
        Vector3d P2_proj =
            PlaneVsSphere::getProjectionOfPoint(P.segment(0, 6), P2);

        double length = (P1 - P2).norm();
        double length_1 = (P1 - P1_proj).norm();
        double length_2 = (P2 - P2_proj).norm();

        double t = 0.5 + curveScale() * length * (length_1 - length_2) /
                             (length_1 + length_2);

        t = sigmoid(t, sigScale());
        Pt_line = P1 + t * (P2 - P1);
        Pt_plane =
            PlaneVsSphere::getProjectionOfPoint(P.segment(0, 6), Pt_line);
    }
    static double sigScale() {
        static double sig_scale = 5.;
        return sig_scale;
    }

    static double curveScale() {
        static double curve_scale = 10.;
        return curve_scale;
    }
};

class PlaneVsPlane {
public:
    static double compute_D(const Vector12d& P, const Vector0d& props) {
        Vector3d normal1 = P.segment(3, 3);
        Vector3d normal2 = P.segment(9, 3);
        if (fabs(1. - fabs(normal1.dot(normal2))) < 1e-8) {
            Vector3d projP =
                PlaneVsSphere::getProjectionOfPoint(P.head(6), P.segment(6, 3));
            return (P.segment(6, 3) - projP).norm();
        } else {
            return 0.;
        }
    }

    FD_CHECK_dDdP(12, 12, 0, 0, "PlaneVsPlane - dDdP_12");
    static void compute_dDdP(Vector12d& dDdP, const Vector12d& P,
                             const Vector0d& props) {
        if (sameNormal(P)) {
            PlaneToPlaneDistance_CodeGen::AD_PlaneToPlaneDistanceGradient(P,
                                                                          dDdP);
        } else {
            dDdP.setZero();
        }
    }
    
    FD_CHECK_d2DdP2(12, 12, 0, 0, "PlaneVsPlane - d2DdP2_12");
    static void compute_d2DdP2(Matrix12d& d2DdP2, const Vector12d& P,
                               const Vector0d& props) {
        if (sameNormal(P)) {
            PlaneToPlaneDistance_CodeGen::AD_PlaneToPlaneDistanceHessian(
                P, d2DdP2);
        } else {
            d2DdP2.setZero();
        }
    }

    FD_CHECK_dDdP(6, 12, 0, 0, "PlaneVsPlane - dDdP_6");
    static void compute_dDdP(Vector6d& dDdP, const Vector12d& P,
                             const Vector0d& props) {
        if (sameNormal(P)) {
            Vector12d dDdP_full;
            PlaneToPlaneDistance_CodeGen::AD_PlaneToPlaneDistanceGradient(
                P, dDdP_full);
            dDdP = dDdP_full.head(6);
        } else {
            dDdP.setZero();
        }
    }

    FD_CHECK_d2DdP2(6, 12, 0, 0, "PlaneVsPlane - d2DdP2_6");
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Vector12d& P,
                               const Vector0d& props) {
        if (sameNormal(P)) {
            Matrix12d d2DdP2_full;
            PlaneToPlaneDistance_CodeGen::AD_PlaneToPlaneDistanceHessian(
                P, d2DdP2_full);
            d2DdP2 = d2DdP2_full.block(0, 0, 6, 6);
        } else {
            d2DdP2.setZero();
        }
    }

private:
    static bool sameNormal(const Vector12d& P) {
        Vector3d n1 = P.segment(3, 3);
        Vector3d n2 = P.segment(9, 3);
        return (fabs(fabs(n1.dot(n2)) - 1.) < 1e-5);
    }
};

}  // namespace DCA

#endif /* __DCA_PLANE_H__ */
#ifndef __DCA__SPHEREVSCAPSULE_H__
#define __DCA__SPHEREVSCAPSULE_H__

#include "AD_PointToLine.h"
#include "FD.h"
#include "utils.h"

namespace DCA {

class SphereVsCapsule {
public:
    static double compute_D(const Vector9d& P, const Vector2d& props) {
        Vector3d P0 = P.head(3);
        return (P0 - computeClosestPointOnLine(P)).norm() - props.sum();
    }

    FD_CHECK_dDdP(9, 9, 2, 0, "SphereVsCapsule - dDdP_9");
    static void compute_dDdP(Vector9d& dDdP, const Vector9d& P,
                             const Vector2d& props) {
        Eigen::Matrix<double, 9, 1> grad;
        PointToLineSegmentDistance_CodeGen::
            AD_PointToLineSegmentDistanceGradient(P, sigScale(), grad);
        dDdP = grad;
    }

    FD_CHECK_dDdP(3, 9, 2, 0, "SphereVsCapsule - dDdP_3");
    static void compute_dDdP(Vector3d& dDdP, const Vector9d& P,
                             const Vector2d& props) {
        Eigen::Matrix<double, 9, 1> grad;
        PointToLineSegmentDistance_CodeGen::
            AD_PointToLineSegmentDistanceGradient(P, sigScale(), grad);

        dDdP = grad.head(3);
    }

    FD_CHECK_dDdP(6, 9, 2, 3, "SphereVsCapsule - dDdP_6");
    static void compute_dDdP(Vector6d& dDdP, const Vector9d& P,
                             const Vector2d& props) {
        Vector9d grad;
        PointToLineSegmentDistance_CodeGen::
            AD_PointToLineSegmentDistanceGradient(P, sigScale(), grad);
        dDdP = grad.tail(6);
    }

    FD_CHECK_d2DdP2(9, 9, 2, 0, "SphereVsCapsule - d2DdP2_9");
    static void compute_d2DdP2(Matrix9d& d2DdP2, const Vector9d& P,
                               const Vector2d& props) {
        Matrix9d hess;
        PointToLineSegmentDistance_CodeGen::
            AD_PointToLineSegmentDistanceHessian(P, sigScale(), hess);
        d2DdP2 = hess;
    }

    FD_CHECK_d2DdP2(3, 9, 2, 0, "SphereVsCapsule - d2DdP2_3");
    static void compute_d2DdP2(Matrix3d& d2DdP2, const Vector9d& P,
                               const Vector2d& props) {
        Matrix9d hess;
        PointToLineSegmentDistance_CodeGen::
            AD_PointToLineSegmentDistanceHessian(P, sigScale(), hess);
        d2DdP2 = hess.block(0, 0, 3, 3);
    }

    FD_CHECK_d2DdP2(6, 9, 2, 3, "SphereVsCapsule - d2DdP2_6");
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Vector9d& P,
                               const Vector2d& props) {
        Matrix9d hess;
        PointToLineSegmentDistance_CodeGen::
            AD_PointToLineSegmentDistanceHessian(P, sigScale(), hess);
        d2DdP2 = hess.block(3, 3, 6, 6);
    }

private:
    static Vector3d computeClosestPointOnLine(const Vector9d& P) {
        Vector3d P1 = P.segment(3, 3);
        Vector3d P2 = P.tail(3);
        return P1 + (P2 - P1) * compute_t(P);
    }

    static double compute_t(const Vector9d& P) {
        Vector3d P0 = P.head(3);
        Vector3d P1 = P.segment(3, 3);
        Vector3d P2 = P.tail(3);
        double t = -1.0 * (((P1 - P0).dot(P2 - P1)) / (P2 - P1).squaredNorm());
        return sigmoid(t, sigScale());
    }

    static double sigScale() {
        static double s = 5.0;
        return s;
    }
};

}  // namespace DCA

#endif
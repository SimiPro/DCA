#include <DCA/Autodiff/AD_PointToLine.h>
#include <DCA/Interactions/SphereVsCapsule.h>

namespace DCA {

double SphereVsCapsule::compute_D(const Vector9d& P, const Vector2d& props) {
    Vector3d P0 = P.head(3);
    return (P0 - computeClosestPointOnLine(P)).norm() - props.sum();
}

void SphereVsCapsule::compute_dDdP(Vector9d& dDdP, const Vector9d& P,
                                   const Vector2d& props) {
    Eigen::Matrix<double, 9, 1> grad;
    PointToLineSegmentDistance_CodeGen::AD_PointToLineSegmentDistanceGradient(
        P, sigScale(), grad);
    dDdP = grad;
}

void SphereVsCapsule::compute_dDdP(Vector3d& dDdP, const Vector9d& P,
                                   const Vector2d& props) {
    Eigen::Matrix<double, 9, 1> grad;
    PointToLineSegmentDistance_CodeGen::AD_PointToLineSegmentDistanceGradient(
        P, sigScale(), grad);

    dDdP = grad.head(3);
}

void SphereVsCapsule::compute_dDdP(Vector6d& dDdP, const Vector9d& P,
                                   const Vector2d& props) {
    Vector9d grad;
    PointToLineSegmentDistance_CodeGen::AD_PointToLineSegmentDistanceGradient(
        P, sigScale(), grad);
    dDdP = grad.tail(6);
}

void SphereVsCapsule::compute_d2DdP2(Matrix9d& d2DdP2, const Vector9d& P,
                                     const Vector2d& props) {
    Matrix9d hess;
    PointToLineSegmentDistance_CodeGen::AD_PointToLineSegmentDistanceHessian(
        P, sigScale(), hess);
    d2DdP2 = hess;
}

void SphereVsCapsule::compute_d2DdP2(Matrix3d& d2DdP2, const Vector9d& P,
                                     const Vector2d& props) {
    Matrix9d hess;
    PointToLineSegmentDistance_CodeGen::AD_PointToLineSegmentDistanceHessian(
        P, sigScale(), hess);
    d2DdP2 = hess.block(0, 0, 3, 3);
}

void SphereVsCapsule::compute_d2DdP2(Matrix6d& d2DdP2, const Vector9d& P,
                                     const Vector2d& props) {
    Matrix9d hess;
    PointToLineSegmentDistance_CodeGen::AD_PointToLineSegmentDistanceHessian(
        P, sigScale(), hess);
    d2DdP2 = hess.block(3, 3, 6, 6);
}

Vector3d SphereVsCapsule::computeClosestPointOnLine(const Vector9d& P) {
    Vector3d P1 = P.segment(3, 3);
    Vector3d P2 = P.tail(3);
    return P1 + (P2 - P1) * compute_t(P);
}

double SphereVsCapsule::compute_t(const Vector9d& P) {
    Vector3d P0 = P.head(3);
    Vector3d P1 = P.segment(3, 3);
    Vector3d P2 = P.tail(3);
    double t = -1.0 * (((P1 - P0).dot(P2 - P1)) / (P2 - P1).squaredNorm());
    return sigmoid(t, sigScale());
}

double SphereVsCapsule::sigScale() {
    static double s = 5.0;
    return s;
}

}  // namespace DCA
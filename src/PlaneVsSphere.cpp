#include <DCA/Autodiff/AD_PlaneToSphere.h>
#include <DCA/Interactions/PlaneVsSphere.h>
#include <DCA/Interactions/SphereVsSphere.h>

namespace DCA {

double PlaneVsSphere::compute_D(const Vector9d& P, const Vector1d& props) {
    Vector2d radius_padded;
    radius_padded << 0, props;

    Vector3d spherePos = P.segment(6, 3);
    Vector3d planePos = getProjectionOfPoint(P.segment(0, 6), spherePos);
    Vector6d bothPos;
    bothPos << planePos, spherePos;
    return SphereVsSphere::compute_D(bothPos, radius_padded);
}

void PlaneVsSphere::compute_dDdP(Vector9d& dDdP, const Vector9d& P,
                                 const Vector1d& props) {
    PlaneToSphereDistance_CodeGen::AD_PlaneToSphereDistanceGradient(P, props,
                                                                    dDdP);
}

void PlaneVsSphere::compute_dDdP(Vector6d& dDdP, const Vector9d& P,
                                 const Vector1d& props) {
    Vector9d dDdP_full;
    PlaneToSphereDistance_CodeGen::AD_PlaneToSphereDistanceGradient(P, props,
                                                                    dDdP_full);
    dDdP = dDdP_full.head(6);
}

void PlaneVsSphere::compute_dDdP(Vector3d& dDdP, const Vector9d& P,
                                 const Vector1d& props) {
    Vector9d dDdP_full;
    PlaneToSphereDistance_CodeGen::AD_PlaneToSphereDistanceGradient(P, props,
                                                                    dDdP_full);
    dDdP = dDdP_full.tail(3);
}

void PlaneVsSphere::compute_d2DdP2(Matrix9d& d2DdP2, const Vector9d& P,
                                   const Vector1d& props) {
    PlaneToSphereDistance_CodeGen::AD_PlaneToSphereDistanceHessian(P, props,
                                                                   d2DdP2);
}

void PlaneVsSphere::compute_d2DdP2(Matrix6d& d2DdP2, const Vector9d& P,
                                   const Vector1d& props) {
    Matrix9d d2DdP2_full;
    PlaneToSphereDistance_CodeGen::AD_PlaneToSphereDistanceHessian(P, props,
                                                                   d2DdP2_full);
    d2DdP2 = d2DdP2_full.block(0, 0, 6, 6);
}

void PlaneVsSphere::compute_d2DdP2(Matrix3d& d2DdP2, const Vector9d& P,
                                   const Vector1d& props) {
    Matrix9d d2DdP2_full;
    PlaneToSphereDistance_CodeGen::AD_PlaneToSphereDistanceHessian(P, props,
                                                                   d2DdP2_full);
    d2DdP2 = d2DdP2_full.block(6, 6, 3, 3);
}

Matrix3d PlaneVsSphere::getPointProjectionDerivative(const Vector9d& P,
                                                     const Vector1d& props) {
    Matrix3d mat = Matrix3d::Identity();
    const Vector3d normal = P.segment(3, 3);
    mat.row(0) += -normal.transpose() * normal(0);
    mat.row(1) += -normal.transpose() * normal(1);
    mat.row(2) += -normal.transpose() * normal(2);
    return mat;
}

Vector3d PlaneVsSphere::getProjectionOfPoint(const Vector6d& P,
                                             const Vector3d& point) {
    return point + -P.segment(3, 3) * getSignedDistanceToPoint(P, point);
}

double PlaneVsSphere::getSignedDistanceToPoint(const Vector6d& P,
                                               const Vector3d& point) {
    return (point - P.head(3)).dot(P.segment(3, 3));
}

void PlaneVsSphere::getCartesianEquationCoeffs(const Vector6d& P, double& a,
                                               double& b, double& c,
                                               double& d) {
    // normal from 3 to 5
    a = P(3);
    b = P(4);
    c = P(5);
    d = -(a * P(0) + b * P(1) + c * P(2));
}
}  // namespace DCA
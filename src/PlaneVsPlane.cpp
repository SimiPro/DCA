#include <DCA/Autodiff/AD_PlaneToPlane.h>
#include <DCA/Interactions/PlaneVsPlane.h>
#include <DCA/Interactions/PlaneVsSphere.h>

namespace DCA {

double PlaneVsPlane::compute_D(const Vector12d& P, const Vector0d& props) {
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

void PlaneVsPlane::compute_dDdP(Vector12d& dDdP, const Vector12d& P,
                                const Vector0d& props) {
    if (sameNormal(P)) {
        PlaneToPlaneDistance_CodeGen::AD_PlaneToPlaneDistanceGradient(P, dDdP);
    } else {
        dDdP.setZero();
    }
}

void PlaneVsPlane::compute_d2DdP2(Matrix12d& d2DdP2, const Vector12d& P,
                                  const Vector0d& props) {
    if (sameNormal(P)) {
        PlaneToPlaneDistance_CodeGen::AD_PlaneToPlaneDistanceHessian(P, d2DdP2);
    } else {
        d2DdP2.setZero();
    }
}

void PlaneVsPlane::compute_dDdP(Vector6d& dDdP, const Vector12d& P,
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

void PlaneVsPlane::compute_d2DdP2(Matrix6d& d2DdP2, const Vector12d& P,
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

bool PlaneVsPlane::sameNormal(const Vector12d& P) {
    Vector3d n1 = P.segment(3, 3);
    Vector3d n2 = P.segment(9, 3);
    return (fabs(fabs(n1.dot(n2)) - 1.) < 1e-5);
}
}  // namespace DCA
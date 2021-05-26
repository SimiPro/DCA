#include <DCA/Autodiff/AD_PointToLine.h>
#include <DCA/Interactions/CapsuleVsCapsule.h>

namespace DCA {

double CapsuleVsCapsule::compute_D(const Vector12d& P, const Vector2d& props) {
    Vector2d X;

    solveForX(P, X);

    return objective().compute_D(P, X) - props.sum();
}

void CapsuleVsCapsule::compute_dDdP(Vector12d& dDdP, const Vector12d& P,
                                    const Vector2d& props) {
    Vector2d X;
    solveForX(P, X);

    Eigen::Matrix<double, 2, 12> dXdP;
    compute_dXdP(dXdP, P, X);

    Vector2d pDpX;
    objective().compute_dDdX(pDpX, P, X);

    Vector12d pDpP;
    objective().compute_dDdP(pDpP, P, X);

    dDdP = dXdP.transpose() * pDpX + pDpP;
}

void CapsuleVsCapsule::compute_dDdP(Vector6d& dDdP, const Vector12d& P,
                                    const Vector2d& props) {
    Vector12d dDdP_full;
    compute_dDdP(dDdP_full, P, props);

    dDdP = dDdP_full.head(6);
}

void CapsuleVsCapsule::compute_d2DdP2(Matrix12d& d2DdP2, const Vector12d& P,
                                      const Vector2d& props) {
    Vector2d X;
    solveForX(P, X);

    Eigen::Matrix<double, 2, 12> dXdP;
    compute_dXdP(dXdP, P, X);

    Matrix2d p2DpX2;
    objective().compute_d2DdX2(p2DpX2, P, X);

    Matrix12d p2DpP2;
    objective().compute_d2DdP2(p2DpP2, P, X);

    Eigen::Matrix<double, 2, 12> p2DpXpP;
    objective().compute_d2DdXdP(p2DpXpP, P, X);

    d2DdP2 = dXdP.transpose() * p2DpX2 * dXdP + dXdP.transpose() * p2DpXpP +
             p2DpXpP.transpose() * dXdP + p2DpP2;
}

void CapsuleVsCapsule::compute_d2DdP2(Matrix6d& d2DdP2, const Vector12d& P,
                                      const Vector2d& props) {
    Matrix12d d2DdP2_full;
    compute_d2DdP2(d2DdP2_full, P, props);

    d2DdP2 = d2DdP2_full.block(0, 0, 6, 6);
}

void CapsuleVsCapsule::compute_dXdP(Eigen::Matrix<double, 2, 12>& dXdP,
                                    const Vector12d& P, const Vector2d& X) {
    Matrix2d dGdX;
    objective().compute_d2OdX2(dGdX, P, X);

    Eigen::Matrix<double, 2, 12> dGdP;
    objective().compute_d2DdXdP(dGdP, P, X);

    dXdP = -1. * dGdX.inverse() * dGdP;
}

void CapsuleVsCapsule::solveForX(const Vector12d& P, Vector2d& X) {
#ifdef RUN_FD_CHECK
    // If we run FD checks, we set the solver residual lower.
    NewtonOptimizer<12, 2> optimizer(1e-12, 15);
#else
    NewtonOptimizer<12, 2> optimizer;

#endif
    X << 0.5, 0.5;
    optimizer.optimize(objective(), P, X, 100);
}

void CapsuleVsCapsule::computeClosestPointOnLines(Vector3d& P12, Vector3d& P34,
                                                  const Vector12d& P) {
    Vector2d X;
    solveForX(P, X);
    Vector3d P1 = P.segment(0, 3);
    Vector3d P2 = P.segment(3, 3);
    Vector3d P3 = P.segment(6, 3);
    Vector3d P4 = P.segment(9, 3);

    P12 = P1 + X[0] * (P2 - P1);
    P34 = P3 + X[1] * (P4 - P3);
}

CapsuleDistanceObjective& CapsuleVsCapsule::objective() {
    static CapsuleDistanceObjective objective;
    return objective;
}

}  // namespace DCA

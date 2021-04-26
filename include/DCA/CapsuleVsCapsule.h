#ifndef __DCA__CAPSULEVSCAPSULE_H__
#define __DCA__CAPSULEVSCAPSULE_H__

#include "AD_PointToLine.h"
#include "CapsuleDistanceObjective.h"
#include "FD.h"
#include "Newton.h"
#include "utils.h"

namespace DCA {

class CapsuleVsCapsule {
public:
    static double compute_D(const Vector12d& P, const Vector2d& props) {
        Vector2d X;

        solveForX(P, X);

        return objective().compute_D(P, X) - props.sum();
    }

    FD_CHECK_dDdP(12, 12, 2, 0, "CapsuleVsCapsule - dDdP_12");
    static void compute_dDdP(Vector12d& dDdP, const Vector12d& P,
                             const Vector2d& props) {
        Vector2d X;
        solveForX(P, X);

        MatrixXd dXdP;
        compute_dXdP(dXdP, P, X);

        VectorXd pDpX;
        objective().compute_dDdX(pDpX, P, X);

        VectorXd pDpP;
        objective().compute_dDdP(pDpP, P, X);

        dDdP = dXdP.transpose() * pDpX + pDpP;
    }

    FD_CHECK_dDdP(6, 12, 2, 0, "CapsuleVsCapsule - dDdP_6");
    static void compute_dDdP(Vector6d& dDdP, const Vector12d& P,
                             const Vector2d& props) {
        Vector12d dDdP_full;
        compute_dDdP(dDdP_full, P, props);

        dDdP = dDdP_full.head(6);
    }

    FD_CHECK_d2DdP2(12, 12, 2, 0, "CapsuleVsCapsule - d2DdP2_12");
    static void compute_d2DdP2(Matrix12d& d2DdP2, const Vector12d& P,
                               const Vector2d& props) {
        Vector2d X;
        solveForX(P, X);

        MatrixXd dXdP;
        compute_dXdP(dXdP, P, X);

        MatrixXd p2DpX2;
        objective().compute_d2DdX2(p2DpX2, P, X);

        MatrixXd p2DpP2;
        objective().compute_d2DdP2(p2DpP2, P, X);

        MatrixXd p2DpXpP;
        objective().compute_d2DdXdP(p2DpXpP, P, X);

        d2DdP2 = dXdP.transpose() * p2DpX2 * dXdP + dXdP.transpose() * p2DpXpP +
                 p2DpXpP.transpose() * dXdP + p2DpP2;
    }

    FD_CHECK_d2DdP2(6, 12, 2, 0, "CapsuleVsCapsule - d2DdP2_6");
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Vector12d& P,
                               const Vector2d& props) {
        Matrix12d d2DdP2_full;
        compute_d2DdP2(d2DdP2_full, P, props);

        d2DdP2 = d2DdP2_full.block(0, 0, 6, 6);
    }

    static void compute_dXdP(MatrixXd& dXdP, const VectorXd& P,
                             const VectorXd& X) {
        dXdP.resize(2, 12);
        MatrixXd d2OdX2;
        objective().compute_d2OdX2(d2OdX2, P, X);

        MatrixXd dGdX = d2OdX2;
        dGdX(0, 1) = dGdX(1, 0);  // is this still needed?

        MatrixXd dGdP;
        objective().compute_d2DdXdP(dGdP, P, X);

        dXdP = -1. * dGdX.inverse() * dGdP;
    }

private:
    static void solveForX(const Vector12d& P, Vector2d& X) {
        NewtonOptimizer optimizer;
        X << 0.5, 0.5;
        // ugly workaround because x here is const size
        VectorXd x_tmp = X;
        optimizer.optimize(objective(), P, x_tmp, 100);
        X = x_tmp;
    }

    /**
     * Compute the closest points on two lines based on P.
     */
    static void computeClosestPointOnLines(Vector3d& P12, Vector3d& P34,
                                           const VectorXd& P) {
        Vector2d X;
        solveForX(P, X);
        Vector3d P1 = P.segment(0, 3);
        Vector3d P2 = P.segment(3, 3);
        Vector3d P3 = P.segment(6, 3);
        Vector3d P4 = P.segment(9, 3);

        P12 = P1 + X[0] * (P2 - P1);
        P34 = P3 + X[1] * (P4 - P3);
    }

    /**
     * Get a reference to the objective
     */
    static CapsuleDistanceObjective& objective() {
        static CapsuleDistanceObjective objective;
        return objective;
    }
};

}  // namespace DCA

#endif
#pragma once

#include <DCA/Autodiff/AD_PlaneToCapsule.h>
#include <DCA/Interactions/PlaneVsSphere.h>
#include <DCA/Utils/FD.h>

namespace DCA {

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

}  // namespace DCA
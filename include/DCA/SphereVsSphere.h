#ifndef __DCA__SPHEREVSSPHERE_H__
#define __DCA__SPHEREVSSPHERE_H__

#include "FD.h"
#include "utils.h"

namespace DCA {

class SphereVsSphere {
public:
    static double compute_D(const Vector6d& P, const Vector2d& props) {
        return (P.head(3) - P.tail(3)).norm() - props.sum();
    }

    FD_CHECK_dDdP(3, 6, 2, 0, "SphereVsSphere - dDdP_3");
    static void compute_dDdP(Vector3d& dDdP, const Vector6d& P,
                             const Vector2d& props) {
        Vector3d v = P.head(3) - P.tail(3);
        double v_norm = v.norm();

        if (v_norm < EPSILON) {
            v_norm = EPSILON;
        }

        dDdP = v / v_norm;
    }

    FD_CHECK_dDdP(6, 6, 2, 0, "SphereVsSphere - dDdP_6");
    static void compute_dDdP(Vector6d& dDdP, const Vector6d& P,
                             const Vector2d& props) {
        Vector3d dDdP1;
        compute_dDdP(dDdP1, P, props);

        dDdP.head(3) = dDdP1;
        dDdP.tail(3) = -dDdP1;
    }

    FD_CHECK_d2DdP2(3, 6, 2, 0, "SphereVsSphere - d2DdP2_3");
    static void compute_d2DdP2(Matrix3d& d2DdP2, const Vector6d& P,
                               const Vector2d& props) {
        Vector3d v = P.head(3) - P.tail(3);
        double v_norm = v.norm();

        if (v_norm < EPSILON) {
            v_norm = EPSILON;
        }
        d2DdP2 =
            (Matrix3d::Identity() * v_norm - (v * v.transpose()) / v_norm) /
            (v_norm * v_norm);
    }

    FD_CHECK_d2DdP2(6, 6, 2, 0, "SphereVsSphere - d2DdP2_6");
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Vector6d& P,
                               const Vector2d& props) {
        Matrix3d d2DdP2_1;
        compute_d2DdP2(d2DdP2_1, P, props);
        d2DdP2.block(0, 0, 3, 3) = d2DdP2_1;
        d2DdP2.block(0, 3, 3, 3) = -d2DdP2_1;
        d2DdP2.block(3, 0, 3, 3) = -d2DdP2_1;
        d2DdP2.block(3, 3, 3, 3) = d2DdP2_1;
    }
};

}  // namespace DCA

#endif
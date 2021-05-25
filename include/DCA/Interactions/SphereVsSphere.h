#ifndef __DCA_SPHEREVSSPHERE_H__
#define __DCA_SPHEREVSSPHERE_H__

#include <DCA/Utils/FD.h>
#include <DCA/Utils/Utils.h>

namespace DCA {

/**
 * @brief This helper is used to compute Sphere vs. Sphere interactions.
 */
class SphereVsSphere {
public:
    /**
     * @brief Compute the distance of two sphere.
     * @param[in] P The parameters for both spheres, that is the degrees of freedom (positions for both spheres), stacked.
     * @param[in] props The properties of both spheres, that is the two radii.
     * @return The distance between both spheres.
     */
    static double compute_D(const Vector6d& P, const Vector2d& props) {
        return (P.head(3) - P.tail(3)).norm() - props.sum();
    }

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    FD_CHECK_dDdP(6, 6, 2, 0, "SphereVsSphere - dDdP_6");
    #endif /* DOXYGEN_SHOULD_SKIP_THIS */
    /**
     * @brief Compute the *full* derivative of the distance of two spheres with respect to P.
     * @param[out] dDdP The *full* derivative \f$\frac{dD}{dP}\f$.
     * @param[in] P The parameters for both spheres, that is the degrees of freedom (positions for both spheres), stacked.
     * @param[in] props The properties of both spheres, that is the two radii.
     */
    static void compute_dDdP(Vector6d& dDdP, const Vector6d& P,
                             const Vector2d& props) {
        Vector3d dDdP1;
        compute_dDdP(dDdP1, P, props);

        dDdP.head(3) = dDdP1;
        dDdP.tail(3) = -dDdP1;
    }

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    FD_CHECK_dDdP(3, 6, 2, 0, "SphereVsSphere - dDdP_3");
    #endif /* DOXYGEN_SHOULD_SKIP_THIS */
    /**
     * @brief Compute the *partial* derivative of the distance of two spheres with respect to P.
     * @param[out] dDdP The *partial* derivative \f$\frac{dD}{dP}\f$, that is with respect to the parameters of the first sphere.
     * @param[in] P The parameters for both spheres, that is the degrees of freedom (positions for both spheres), stacked.
     * @param[in] props The properties of both spheres, that is the two radii.
     */
    static void compute_dDdP(Vector3d& dDdP, const Vector6d& P,
                             const Vector2d& props) {
        Vector3d v = P.head(3) - P.tail(3);
        double v_norm = v.norm();

        if (v_norm < EPSILON()) {
            v_norm = EPSILON();
        }

        dDdP = v / v_norm;
    }

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    FD_CHECK_d2DdP2(6, 6, 2, 0, "SphereVsSphere - d2DdP2_6");
    #endif /* DOXYGEN_SHOULD_SKIP_THIS */
    /**
     * @brief Compute the *full* second derivative of the distance of two spheres with respect to P.
     * @param[out] d2DdP2 The *full* second derivative \f$\frac{d^2D}{dP^2}\f$.
     * @param[in] P The parameters for both spheres, that is the degrees of freedom (positions for both spheres), stacked.
     * @param[in] props The properties of both spheres, that is the two radii.
     */
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Vector6d& P,
                               const Vector2d& props) {
        Matrix3d d2DdP2_1;
        compute_d2DdP2(d2DdP2_1, P, props);
        d2DdP2.block(0, 0, 3, 3) = d2DdP2_1;
        d2DdP2.block(0, 3, 3, 3) = -d2DdP2_1;
        d2DdP2.block(3, 0, 3, 3) = -d2DdP2_1;
        d2DdP2.block(3, 3, 3, 3) = d2DdP2_1;
    }

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    FD_CHECK_d2DdP2(3, 6, 2, 0, "SphereVsSphere - d2DdP2_3");
    #endif /* DOXYGEN_SHOULD_SKIP_THIS */
    /**
     * @brief Compute the *partial* second derivative of the distance of two spheres with respect to P.
     * @param[out] d2DdP2 The *partial* derivative \f$\frac{d^2D}{dP^2}\f$, that is with respect to the parameters of the first sphere.
     * @param[in] P The parameters for both spheres, that is the degrees of freedom (positions for both spheres), stacked.
     * @param[in] props The properties of both spheres, that is the two radii.
     */
    static void compute_d2DdP2(Matrix3d& d2DdP2, const Vector6d& P,
                               const Vector2d& props) {
        Vector3d v = P.head(3) - P.tail(3);
        double v_norm = v.norm();

        if (v_norm < EPSILON()) {
            v_norm = EPSILON();
        }
        d2DdP2 =
            (Matrix3d::Identity() * v_norm - (v * v.transpose()) / v_norm) /
            (v_norm * v_norm);
    }
private:
    /**
     * @brief Returns a small number.
     * 
     * Used internally.
     * @return A small number.
     */
    static double EPSILON() {
        return 1e-8;
    }
};

}  // namespace DCA

#endif /* __DCA_SPHEREVSSPHERE_H__ */
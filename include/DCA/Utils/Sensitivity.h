#pragma once

#include <DCA/Utils/Utils.h>

namespace DCA {

/**
 * @brief This base class holds all functions needed in sensitivity analysis.
 * @param SizeP Size of the parameter vector.
 * @param SizeX Size of the solving variable vector.
 */
template <int SizeP, int SizeX>
class SensitivityObjective {
public:
    using P_v = Eigen::Matrix<double, SizeP, 1>;      ///< Helper: P Vector
    using X_v = Eigen::Matrix<double, SizeX, 1>;      ///< Helper: X Vector
    using X_m = Eigen::Matrix<double, SizeX, SizeX>;  ///< Helper: X Matrix
    using P_m = Eigen::Matrix<double, SizeP, SizeP>;  ///< Helper: P Matrix
    using XP_m =
        Eigen::Matrix<double, SizeX, SizeP>;  ///< Helper: X by P Matrix
    /**
     * @brief Computes the distance for a given P and X.
     * @param[in] P Some parameters for the objective.
     * @param[in] X The solving variable.
     * @return The distance.
     */
    virtual double compute_D(const P_v& P, const X_v& X) const = 0;

    /**
     * @brief Computes the first derivative of the distance with respect to X.
     * @param[out] dDdX The derivative \f$\frac{dD}{dX}\f$ of the distance.
     * @param[in] P Some parameters for the objective.
     * @param[in] X The solving variable.
     */
    virtual void compute_dDdX(X_v& dDdX, const P_v& P, const X_v& X) const = 0;

    /**
     * @brief Computes the second derivative of the distance with respect to X.
     * @param[out] d2DdX2 The second derivative \f$\frac{d^2D}{dX^2}\f$ of the distance.
     * @param[in] P Some parameters for the objective.
     * @param[in] X The solving variable.
     */
    virtual void compute_d2DdX2(X_m& d2DdX2, const P_v& P,
                                const X_v& X) const = 0;

    /**
     * @brief Computes the first derivative of the distance with respect to P.
     * @param[out] dDdP The second derivative \f$\frac{dD}{dP}\f$ of the distance.
     * @param[in] P Some parameters for the objective.
     * @param[in] X The solving variable.
     */
    virtual void compute_dDdP(P_v& dDdP, const P_v& P, const X_v& X) const = 0;
    /**
     * @brief Computes the second derivative of the distance with respect to P.
     * @param[out] d2DdP2 The second derivative \f$\frac{d^2D}{dP^2}\f$ of the distance.
     * @param[in] P Some parameters for the objective.
     * @param[in] X The solving variable.
     */
    virtual void compute_d2DdP2(P_m& d2DdP2, const P_v& P,
                                const X_v& X) const = 0;
    /**
     * @brief Computes the second derivative of the distance with respect to X and P.
     * @param[out] d2DdXdP The mixed second derivative \f$\frac{d^2D}{dXdP}\f$ of the distance.
     * @param[in] P Some parameters for the objective.
     * @param[in] X The solving variable.
     */
    virtual void compute_d2DdXdP(XP_m& d2DdXdP, const P_v& P,
                                 const X_v& X) const = 0;
};

}  // namespace DCA
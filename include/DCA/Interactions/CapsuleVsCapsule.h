#pragma once

#include <DCA/Interactions/CapsuleDistanceObjective.h>
#include <DCA/Utils/FD.h>
#include <DCA/Utils/Newton.h>

namespace DCA {

/**
 * @brief This helper is used to compute Capsule vs. Capsule interactions.
 */
class CapsuleVsCapsule {
public:
    /**
     * @brief Compute the distance of two capsules.
     * @param[in] P The parameters for both capsules, that is the degrees of freedom (start- and end-positions for both capsules), stacked.
     * @param[in] props The properties of both capsules, that is the two radii.
     * @return The distance between both capsules.
     */
    static double compute_D(const Vector12d& P, const Vector2d& props);

    FD_CHECK_dDdP(12, 12, 2, 0, "CapsuleVsCapsule - dDdP_12");
    /**
     * @brief Compute the *full* derivative of the distance of two capsules with respect to P.
     * @param[out] dDdP The *full* derivative \f$\frac{dD}{dP}\f$.
     * @param[in] P The parameters for both capsules, that is the degrees of freedom (start- and end-positions for both capsules), stacked.
     * @param[in] props The properties of both capsules, that is the two radii.
     */
    static void compute_dDdP(Vector12d& dDdP, const Vector12d& P,
                             const Vector2d& props);

    FD_CHECK_dDdP(6, 12, 2, 0, "CapsuleVsCapsule - dDdP_6");
    /**
     * @brief Compute the *partial* derivative of the distance of two capsules with respect to P.
     * @param[out] dDdP The *partial* derivative \f$\frac{dD}{dP}\f$, that is with respect to the parameters of the first capsule.
     * @param[in] P The parameters for both capsules, that is the degrees of freedom (start- and end-positions for both capsules), stacked.
     * @param[in] props The properties of both capsules, that is the two radii.
     */
    static void compute_dDdP(Vector6d& dDdP, const Vector12d& P,
                             const Vector2d& props);

    FD_CHECK_d2DdP2(12, 12, 2, 0, "CapsuleVsCapsule - d2DdP2_12");
    /**
     * @brief Compute the *full* second derivative of the distance of two capsules with respect to P.
     * @param[out] d2DdP2 The *full* second derivative \f$\frac{d^2D}{dP^2}\f$.
     * @param[in] P The parameters for both capsules, that is the degrees of freedom (start- and end-positions for both capsules), stacked.
     * @param[in] props The properties of both capsules, that is the two radii.
     */
    static void compute_d2DdP2(Matrix12d& d2DdP2, const Vector12d& P,
                               const Vector2d& props);

    FD_CHECK_d2DdP2(6, 12, 2, 0, "CapsuleVsCapsule - d2DdP2_6");
    /**
     * @brief Compute the *partial* second derivative of the distance of two capsules with respect to P.
     * @param[out] d2DdP2 The *partial* second derivative \f$\frac{d^2D}{dP^2}\f$, that is with respect to the parameters of the first capsule.
     * @param[in] P The parameters for both capsules, that is the degrees of freedom (start- and end-positions for both capsules), stacked.
     * @param[in] props The properties of both capsules, that is the two radii.
     */
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Vector12d& P,
                               const Vector2d& props);

    FD_CHECK_dXdP(12, 12, 2, 0, "CapsuleVsCapsule - dXdP_12");
    /**
     * @brief Compute the derivative of the current point on both lines with respect to P.
     * @param[out] dXdP The *partial* derivative \f$\frac{dX}{dP}\f$, that is with respect to the parameters of both capsule.
     * @param[in] P The parameters for both capsules, that is the degrees of freedom (start- and end-positions for both capsules), stacked.
     * @param[out] X The solution, which are the two points on both lines.
     */
    static void compute_dXdP(Eigen::Matrix<double, 2, 12>& dXdP,
                             const Vector12d& P, const Vector2d& X);

#if DCA_DEV_TOOLS == 1
public:
#else
private:
#endif

    /**
     * @brief Solve the optimization function for X given P.
     * @param[in] P The parameters for both capsules.
     * @param[out] X The solution, which are the two points on both lines.
     */
    static void solveForX(const Vector12d& P, Vector2d& X);

    /**
     * @brief Compute the closest points on two lines based on P.
     * @param[out] P12 The point on the line 1-2
     * @param[out] P34 The point on the line 3-4
     * @param[in] P The parameters for both capsules.
     * @note Relies on solveForX
     */
    static void computeClosestPointOnLines(Vector3d& P12, Vector3d& P34,
                                           const Vector12d& P);

    /**
     * @brief Get a reference to the objective
     * @return A reference to the distance objective.
     */
    static CapsuleDistanceObjective& objective();
};

}  // namespace DCA
#pragma once

#include <DCA/Utils/FD.h>

namespace DCA {

/**
 * @brief This helper is used to compute Plane vs. Plane interactions.
 */
class PlaneVsPlane {
public:
    /**
     * @brief Compute the distance between two planes.
     * @param[in] P The parameters for both planes, that is the degrees of freedom (position and normal for each plane), stacked.
     * @param[in] props The properties of both primitives, that is nothing.
     * @return The distance between the two planes.
     */
    static double compute_D(const Vector12d& P, const Vector0d& props);

    FD_CHECK_dDdP(12, 12, 0, 0, "PlaneVsPlane - dDdP_12");
    /**
     * @brief Compute the *full* derivative of the distance between two planes with respect to P.
     * @param[out] dDdP The *full* derivative \f$\frac{dD}{dP}\f$.
     * @param[in] P The parameters for both planes, that is the degrees of freedom (position and normal for each plane), stacked.
     * @param[in] props The properties of both primitives, that is nothing.
     */
    static void compute_dDdP(Vector12d& dDdP, const Vector12d& P,
                             const Vector0d& props);

    FD_CHECK_d2DdP2(12, 12, 0, 0, "PlaneVsPlane - d2DdP2_12");
    /**
     * @brief Compute the *full* derivative of the distance between two planes with respect to P.
     * @param[out] dDdP The *full* derivative \f$\frac{dD}{dP}\f$.
     * @param[in] P The parameters for both planes, that is the degrees of freedom (position and normal for each plane), stacked.
     * @param[in] props The properties of both primitives, that is nothing.
     */
    static void compute_d2DdP2(Matrix12d& d2DdP2, const Vector12d& P,
                               const Vector0d& props);

    FD_CHECK_dDdP(6, 12, 0, 0, "PlaneVsPlane - dDdP_6");
    /**
     * @brief Compute the *partial* derivative of the distance between two planes with respect to P.
     * @param[out] dDdP The *partial* derivative \f$\frac{dD}{dP}\f$, that is with respect to the parameters of the first plane.
     * @param[in] P The parameters for both planes, that is the degrees of freedom (position and normal for each plane), stacked.
     * @param[in] props The properties of both primitives, that is nothing.
     */
    static void compute_dDdP(Vector6d& dDdP, const Vector12d& P,
                             const Vector0d& props);

    FD_CHECK_d2DdP2(6, 12, 0, 0, "PlaneVsPlane - d2DdP2_6");
    /**
     * @brief Compute the *partial* second derivative of the distance between two planes with respect to P.
     * @param[out] d2DdP2 The *partial* second derivative \f$\frac{d^2D}{dP^2}\f$, that is with respect to the parameters of the sphere.
     * @param[in] P The parameters for both planes, that is the degrees of freedom (position and normal for each plane), stacked.
     * @param[in] props The properties of both primitives, that is nothing.
     */
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Vector12d& P,
                               const Vector0d& props);

private:
    /**
     * @brief Returns whether two planes have the same normal.
     * @param[in] P The parameters for both planes, that is the degrees of freedom (position and normal for each plane), stacked.
     * @return True if both planes have the same normal (also if the normals are flipped).
     */
    static bool sameNormal(const Vector12d& P);
};

}  // namespace DCA

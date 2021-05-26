#pragma once

#include <DCA/Utils/FD.h>

namespace DCA {

/**
 * @brief This helper is used to compute Sphere vs. Capsule interactions.
 */
class SphereVsCapsule {
public:
    /**
     * @brief Compute the distance between a sphere and a capsule.
     * @param[in] P The parameters for the sphere and the capsule, that is the degrees of freedom (position of the sphere, start- and end-positions for the capsule), stacked.
     * @param[in] props The properties of both primitives, that is the two radii.
     * @return The distance between the sphere and the capsule.
     */
    static double compute_D(const Vector9d& P, const Vector2d& props);

    FD_CHECK_dDdP(9, 9, 2, 0, "SphereVsCapsule - dDdP_9");
    /**
     * @brief Compute the *full* derivative of the distance between a sphere and a capsule with respect to P.
     * @param[out] dDdP The *full* derivative \f$\frac{dD}{dP}\f$.
     * @param[in] P The parameters for the sphere and the capsule, that is the degrees of freedom (position of the sphere, start- and end-positions for the capsule), stacked.
     * @param[in] props The properties of both primitives, that is the two radii.
     */
    static void compute_dDdP(Vector9d& dDdP, const Vector9d& P,
                             const Vector2d& props);

    FD_CHECK_dDdP(3, 9, 2, 0, "SphereVsCapsule - dDdP_3");
    /**
     * @brief Compute the *partial* derivative of the distance between a sphere and a capsule with respect to P.
     * @param[out] dDdP The *partial* derivative \f$\frac{dD}{dP}\f$, that is with respect to the parameters of the sphere.
     * @param[in] P The parameters for the sphere and the capsule, that is the degrees of freedom (position of the sphere, start- and end-positions for the capsule), stacked.
     * @param[in] props The properties of both primitives, that is the two radii.
     */
    static void compute_dDdP(Vector3d& dDdP, const Vector9d& P,
                             const Vector2d& props);

    FD_CHECK_dDdP(6, 9, 2, 3, "SphereVsCapsule - dDdP_6");
    /**
     * @brief Compute the *partial* derivative of the distance between a sphere and a capsule with respect to P.
     * @param[out] dDdP The *partial* derivative \f$\frac{dD}{dP}\f$, that is with respect to the parameters of the capsule.
     * @param[in] P The parameters for the sphere and the capsule, that is the degrees of freedom (position of the sphere, start- and end-positions for the capsule), stacked.
     * @param[in] props The properties of both primitives, that is the two radii.
     */
    static void compute_dDdP(Vector6d& dDdP, const Vector9d& P,
                             const Vector2d& props);

    FD_CHECK_d2DdP2(9, 9, 2, 0, "SphereVsCapsule - d2DdP2_9");
    /**
     * @brief Compute the *full* second derivative of the distance between a sphere and a capsule with respect to P.
     * @param[out] d2DdP2 The *full* second derivative \f$\frac{d^2D}{dP^2}\f$.
     * @param[in] P The parameters for the sphere and the capsule, that is the degrees of freedom (position of the sphere, start- and end-positions for the capsule), stacked.
     * @param[in] props The properties of both primitives, that is the two radii.
     */
    static void compute_d2DdP2(Matrix9d& d2DdP2, const Vector9d& P,
                               const Vector2d& props);

    FD_CHECK_d2DdP2(3, 9, 2, 0, "SphereVsCapsule - d2DdP2_3");
    /**
     * @brief Compute the *partial* second derivative of the distance between a sphere and a capsule with respect to P.
     * @param[out] d2DdP2 The *partial* second derivative \f$\frac{d^2D}{dP^2}\f$, that is with respect to the parameters of the sphere.
     * @param[in] P The parameters for the sphere and the capsule, that is the degrees of freedom (position of the sphere, start- and end-positions for the capsule), stacked.
     * @param[in] props The properties of both primitives, that is the two radii.
     */
    static void compute_d2DdP2(Matrix3d& d2DdP2, const Vector9d& P,
                               const Vector2d& props);

    FD_CHECK_d2DdP2(6, 9, 2, 3, "SphereVsCapsule - d2DdP2_6");
    /**
     * @brief Compute the *partial* second derivative of the distance between a sphere and a capsule with respect to P.
     * @param[out] d2DdP2 The *partial* second derivative \f$\frac{d^2D}{dP^2}\f$, that is with respect to the parameters of the capsule.
     * @param[in] P The parameters for the sphere and the capsule, that is the degrees of freedom (position of the sphere, start- and end-positions for the capsule), stacked.
     * @param[in] props The properties of both primitives, that is the two radii.
     */
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Vector9d& P,
                               const Vector2d& props);

private:
    /**
     * @brief Compute the closest point on a line (distance to a point)
     * @param[in] P Stacked representation of the external point (first three entries of P) and the line (start- and end-position, last 6 entries of P).
     * @return The closest point on the line to the external point.
     */
    static Vector3d computeClosestPointOnLine(const Vector9d& P);

    /**
     * @brief Compute the value on the line given by the last six entries in P
     * 
     * The value on the line is the percentage on the line which is closest to the external point (the first three entries in P).
     * @param[in] P Stacked representation of the external point (first three entries of P) and the line (start- and end-position, last 6 entries of P).
     * @return The percentage on the line.
     */
    static double compute_t(const Vector9d& P);
    /**
     * @brief Helper to push values to the center.
     * @return A constant, 5.
     */
    static double sigScale();
};

}  // namespace DCA

#pragma once

#include <DCA/Utils/FD.h>
#include <DCA/Utils/Utils.h>

namespace DCA {

/**
 * @brief This helper is used to compute Plane vs. Capsule interactions.
 */
class PlaneVsCapsule {
public:
    /**
     * @brief Compute the distance between a plane and a capsule.
     * @param[in] P The parameters for the plane and the capsule, that is the degrees of freedom (position and normal of the plane, start- and end-positions for the capsule), stacked.
     * @param[in] props The properties of both primitives, that is the radius of the capsule.
     * @return The distance between the plane and the capsule.
     */
    static double compute_D(const Vector12d& P, const Vector1d& props);

    FD_CHECK_dDdP(6, 12, 1, 0, "PlaneVsCapsule - dDdP_6");
    /**
     * @brief Compute the *partial* derivative of the distance between a plane and a capsule with respect to P.
     * @param[out] dDdP The *partial* derivative \f$\frac{dD}{dP}\f$, that is with respect to the parameters of the plane.
     * @param[in] P The parameters for the plane and the capsule, that is the degrees of freedom (position and normal of the plane, start- and end-positions for the capsule), stacked.
     * @param[in] props The properties of both primitives, that is the radius of the capsule.
     */
    static void compute_dDdP(Vector6d& dDdP, const Vector12d& P,
                             const Vector1d& props);

    FD_CHECK_dDdP(12, 12, 1, 0, "PlaneVsCapsule - dDdP_12");
    /**
     * @brief Compute the *full* derivative of the distance between a plane and a capsule with respect to P.
     * @param[out] dDdP The *full* derivative \f$\frac{dD}{dP}\f$.
     * @param[in] P The parameters for the plane and the capsule, that is the degrees of freedom (position and normal of the plane, start- and end-positions for the capsule), stacked.
     * @param[in] props The properties of both primitives, that is the radius of the capsule.
     */
    static void compute_dDdP(Vector12d& dDdP, const Vector12d& P,
                             const Vector1d& props);

    FD_CHECK_d2DdP2(12, 12, 1, 0, "PlaneVsCapsule - d2DdP2_12");
    /**
     * @brief Compute the *full* second derivative of the distance between a plane and a capsule with respect to P.
     * @param[out] d2DdP2 The *full* second derivative \f$\frac{d^2D}{dP^2}\f$.
     * @param[in] P The parameters for the plane and the capsule, that is the degrees of freedom (position and normal of the plane, start- and end-positions for the capsule), stacked.
     * @param[in] props The properties of both primitives, that is the radius of the capsule.
     */
    static void compute_d2DdP2(Matrix12d& d2DdP2, const Vector12d& P,
                               const Vector1d& props);

    FD_CHECK_d2DdP2(6, 12, 1, 0, "PlaneVsCapsule - d2DdP2_6");
    /**
     * @brief Compute the *partial* second derivative of the distance between a plane and a capsule with respect to P.
     * @param[out] d2DdP2 The *partial* second derivative \f$\frac{d^2D}{dP^2}\f$, that is with respect to the parameters of the plane.
     * @param[in] P The parameters for the plane and the capsule, that is the degrees of freedom (position and normal of the plane, start- and end-positions for the capsule), stacked.
     * @param[in] props The properties of both primitives, that is the radius of the capsule.
     */
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Vector12d& P,
                               const Vector1d& props);

public:
    /**
     * @brief Compute the closest point on the plane and on the capsule.
     * @param[in] P The parameters for the plane and the capsule, that is the degrees of freedom (position and normal of the plane, start- and end-positions for the capsule), stacked.
     * @param[out] Pt_plane 
     * @param[out] Pt_line 
    */
    static void computeClosestPoints(const Vector12d& P, Vector3d& Pt_plane,
                                     Vector3d& Pt_line);
    /**
     * @brief Helper to push values to the center.
     * @return A constant, 5.
     */
    static double sigScale();
    /**
     * @brief Helper to push values to the center.
     * @return A constant, 10.
     */
    static double curveScale();
};

}  // namespace DCA
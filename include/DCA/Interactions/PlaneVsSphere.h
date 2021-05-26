#pragma once

#include <DCA/Utils/FD.h>

namespace DCA {

/**
 * @brief This helper is used to compute Plane vs. Sphere interactions.
 */
class PlaneVsSphere {
public:
    /**
     * @brief Compute the distance between a plane and a sphere.
     * @param[in] P The parameters for the plane and the sphere, that is the degrees of freedom (position and normal of the plane, position of the sphere), stacked.
     * @param[in] props The properties of both primitives, that is the radius of the sphere.
     * @return The distance between the plane and the sphere.
     */
    static double compute_D(const Vector9d& P, const Vector1d& props);

    FD_CHECK_dDdP(9, 9, 1, 0, "PlaneVsSphere - dDdP_9");
    /**
     * @brief Compute the *full* derivative of the distance between a plane and a sphere with respect to P.
     * @param[out] dDdP The *full* derivative \f$\frac{dD}{dP}\f$.
     * @param[in] P The parameters for the plane and the sphere, that is the degrees of freedom (position and normal of the plane, position of the sphere), stacked.
     * @param[in] props The properties of both primitives, that is the radius of the sphere.
     */
    static void compute_dDdP(Vector9d& dDdP, const Vector9d& P,
                             const Vector1d& props);

    FD_CHECK_dDdP(6, 9, 1, 0, "PlaneVsSphere - dDdP_6");
    /**
     * @brief Compute the *partial* derivative of the distance between a plane and a sphere with respect to P.
     * @param[out] dDdP The *partial* derivative \f$\frac{dD}{dP}\f$, that is with respect to the parameters of the plane.
     * @param[in] P The parameters for the plane and the sphere, that is the degrees of freedom (position and normal of the plane, position of the sphere), stacked.
     * @param[in] props The properties of both primitives, that is the radius of the sphere.
     */
    static void compute_dDdP(Vector6d& dDdP, const Vector9d& P,
                             const Vector1d& props);

    FD_CHECK_dDdP(3, 9, 1, 6, "PlaneVsSphere - dDdP_3");
    /**
     * @brief Compute the *partial* derivative of the distance between a plane and a sphere with respect to P.
     * @param[out] dDdP The *partial* derivative \f$\frac{dD}{dP}\f$, that is with respect to the parameters of the sphere.
     * @param[in] P The parameters for the plane and the sphere, that is the degrees of freedom (position and normal of the plane, position of the sphere), stacked.
     * @param[in] props The properties of both primitives, that is the radius of the sphere.
     */
    static void compute_dDdP(Vector3d& dDdP, const Vector9d& P,
                             const Vector1d& props);

    FD_CHECK_d2DdP2(9, 9, 1, 0, "PlaneVsSphere - d2DdP2_9");
    /**
     * @brief Compute the *full* second derivative of the distance between a plane and a sphere with respect to P.
     * @param[out] d2DdP2 The *full* second derivative \f$\frac{d^2D}{dP^2}\f$.
     * @param[in] P The parameters for the plane and the sphere, that is the degrees of freedom (position and normal of the plane, position of the sphere), stacked.
     * @param[in] props The properties of both primitives, that is the radius of the sphere.
     */
    static void compute_d2DdP2(Matrix9d& d2DdP2, const Vector9d& P,
                               const Vector1d& props);
    FD_CHECK_d2DdP2(6, 9, 1, 0, "PlaneVsSphere - d2DdP2_6");
    /**
     * @brief Compute the *partial* second derivative of the distance between a plane and a sphere with respect to P.
     * @param[out] d2DdP2 The *partial* second derivative \f$\frac{d^2D}{dP^2}\f$, that is with respect to the parameters of the plane.
     * @param[in] P The parameters for the plane and the sphere, that is the degrees of freedom (position and normal of the plane, position of the sphere), stacked.
     * @param[in] props The properties of both primitives, that is the radius of the sphere.
     */
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Vector9d& P,
                               const Vector1d& props);
    FD_CHECK_d2DdP2(3, 9, 1, 6, "PlaneVsSphere - d2DdP2_3");
    /**
     * @brief Compute the *partial* second derivative of the distance between a plane and a sphere with respect to P.
     * @param[out] d2DdP2 The *partial* second derivative \f$\frac{d^2D}{dP^2}\f$, that is with respect to the parameters of the sphere.
     * @param[in] P The parameters for the plane and the sphere, that is the degrees of freedom (position and normal of the plane, position of the sphere), stacked.
     * @param[in] props The properties of both primitives, that is the radius of the sphere.
     */
    static void compute_d2DdP2(Matrix3d& d2DdP2, const Vector9d& P,
                               const Vector1d& props);

public:
    /**
     * @brief Compute the derivative of getProjectionOfPoint.
     * @param[in] P The parameters for the plane and the sphere, that is the degrees of freedom (position and normal of the plane, position of the sphere), stacked.
     * @param[in] props The properties of both primitives, that is the radius of the sphere.
     * @return The full derivative of the projection of the point.
     */
    static Matrix3d getPointProjectionDerivative(const Vector9d& P,
                                                 const Vector1d& props);

    /**
     * @brief Get a projection of a point onto a plane.
     * @param[in] P The parameters for the plane, that is the degrees of freedom (position and normal of the plane), stacked.
     * @param[in] point The point to project onto the plane.
     * @return The projected point.
     */
    static Vector3d getProjectionOfPoint(const Vector6d& P,
                                         const Vector3d& point);

    /**
     * @brief Get a signed distance for a point to the plane.
     * @param[in] P The parameters for the plane, that is the degrees of freedom (position and normal of the plane), stacked.
     * @param[in] point The point to get the distance to the plane.
     * @return The signed distance from the plane to the point.
     */
    static double getSignedDistanceToPoint(const Vector6d& P,
                                           const Vector3d& point);

    /**
     * @brief Get the coefficients for the plane equation, that is \f$ \alpha x+\beta y+\gamma c=\delta\f$.
     * @param[in] P The parameters for the plane, that is the degrees of freedom (position and normal of the plane), stacked.
     * @param[out] a The coefficient alpha.
     * @param[out] b The coefficient beta.
     * @param[out] c The coefficient gamma.
     * @param[out] d The coefficient delta.
     */
    static void getCartesianEquationCoeffs(const Vector6d& P, double& a,
                                           double& b, double& c, double& d);
};

}  // namespace DCA
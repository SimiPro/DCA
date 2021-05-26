#pragma once

#include <DCA/Pair.h>
#include <DCA/Utils/Primitives.h>

namespace DCA {

/**
 * @brief Public %API.
 * 
 * This is the main public %API which should be used.
 */
class API {
public:
    /**
     * @defgroup SphereVsSphere
     * @brief %API for %Sphere vs. %Sphere interactions.
     * 
     * This module contains all helpers for %Sphere vs. %Sphere interactions.
     * @{
     */

    /** 
     * @brief Computes the distance between two spheres.
     * 
     * @param s_a The first sphere
     * @param s_b The second sphere
     * @return The distance between both spheres.
     */
    static double compute_D(const Sphere& s_a, const Sphere& s_b);

    /**
     * @brief Computes the *full* derivative of the distance between both spheres.
     * 
     * The derivative is of size 6,
     * meaning it is taken with respect to the parameters of both spheres (both positions, a total of 6 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for both spheres.
     * @param s_a The first sphere
     * @param s_b The second sphere
     */
    static void compute_dDdP(Vector6d& dDdP, const Sphere& s_a,
                             const Sphere& s_b);

    /**
     * @brief Computes the *partial* derivative of the distance between both spheres.
     * 
     * The derivative is of size 3,
     * meaning it is taken with respect to the parameters of the first sphere (namely the position, a total of 3 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for the first sphere.
     * @param s_a The first sphere
     * @param s_b The second sphere
     * 
     * @attention This is not the *full* derivative, only the segment which belong to the first sphere.
     */
    static void compute_dDdP(Vector3d& dDdP, const Sphere& s_a,
                             const Sphere& s_b);

    /**
     * @brief Computes the *full* second derivative of the distance between both spheres.
     * 
     * The derivative is of size 6 by 6,
     * meaning it is taken with respect to the parameters of both spheres (both positions, a total of 6 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for both spheres.
     * @param s_a The first sphere
     * @param s_b The second sphere
     */
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Sphere& s_a,
                               const Sphere& s_b);

    /**
     * @brief Computes the *partial* second derivative of the distance between both spheres.
     * 
     * The derivative is of size 3 by 3,
     * meaning it is taken with respect to the parameters of the first sphere (namely the position, a total of 3 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for the first sphere.
     * @param s_a The first sphere
     * @param s_b The second sphere
     * 
     * @attention This is not the *full* derivative, only the block belonging to the first sphere.
     */
    static void compute_d2DdP2(Matrix3d& d2DdP2, const Sphere& s_a,
                               const Sphere& s_b);

    /**
     * @}
     */ // end of group SphereVsSphere

    /**
     * @defgroup SphereVsCapsule
     * @brief %API for %Sphere vs. %Capsule interactions.
     * 
     * This module contains all helpers for %Sphere vs. %Capsule interactions.
     * @{
     */

    /** 
     * @brief Computes the distance between a sphere and a capsule.
     * 
     * @param s_a The sphere
     * @param c_b The capsule
     * @return The distance between both primitives.
     */
    static double compute_D(const Sphere& s_a, const Capsule& c_b);

    /**
     * @brief Computes the *full* derivative of the distance between a sphere and a capsule.
     * 
     * The derivative is size 9,
     * meaning it is taken with respect to the parameters of both primitives (both positions, a total of 9 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for both primitives.
     * @param s_a The sphere
     * @param c_b The capsule
     */
    static void compute_dDdP(Vector9d& dDdP, const Sphere& s_a,
                             const Capsule& c_b);

    /**
     * @brief Computes the *partial* derivative of the distance between a sphere and a capsule.
     * 
     * The derivative is of size 3,
     * meaning it is taken with respect to the parameters of the sphere (position, a total of 3 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for the sphere.
     * @param s_a The sphere
     * @param c_b The capsule
     * 
     * @attention This is not the *full* derivative, only the segment which belongs to the sphere.
     */
    static void compute_dDdP(Vector3d& dDdP, const Sphere& s_a,
                             const Capsule& c_b);

    /**
     * @brief Computes the *partial* derivative of the distance between a sphere and a capsule.
     * 
     * The derivative is of size 6,
     * meaning it is taken with respect to the parameters of the capsule (start- and end-position, a total of 6 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for the capsule.
     * @param s_a The sphere
     * @param c_b The capsule
     * 
     * @attention This is not the *full* derivative, only the segment which belongs to the capsule.
     */
    static void compute_dDdP(Vector6d& dDdP, const Sphere& s_a,
                             const Capsule& c_b);

    /**
     * @brief Computes the *full* second derivative of the distance between a sphere and a capsule.
     * 
     * The derivative is size 9 by 9,
     * meaning it is taken with respect to the parameters of both primitives (all positions, a total of 9 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP2} \f$ where D is the distance and P is the vector of parameters for both primitives.
     * @param s_a The sphere
     * @param c_b The capsule
     */
    static void compute_d2DdP2(Matrix9d& d2DdP2, const Sphere& s_a,
                               const Capsule& c_b);

    /**
     * @brief Computes the *partial* second derivative of the distance between a sphere and a capsule.
     * 
     * The derivative is of size 3 by 3,
     * meaning it is taken with respect to the parameters of the sphere (position, a total of 3 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for the sphere.
     * @param s_a The sphere
     * @param c_b The capsule
     * 
     * @attention This is not the *full* derivative, only the block which belongs to the sphere.
     */
    static void compute_d2DdP2(Matrix3d& d2DdP2, const Sphere& s_a,
                               const Capsule& c_b);

    /**
     * @brief Computes the *partial* second derivative of the distance between a sphere and a capsule.
     * 
     * The derivative is of size 6 by 6,
     * meaning it is taken with respect to the parameters of the capsule (start- and end-position, a total of 6 degrees of freedom).
     * 
     * @param[out] d2DdP2 The first derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for the capsule.
     * @param s_a The sphere
     * @param c_b The capsule
     * 
     * @attention This is not the *full* second derivative, only the block which belongs to the capsule.
     */
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Sphere& s_a,
                               const Capsule& c_b);

    /**
     * @}
     */ // end of group SphereVsCapsule

    /**
     * @defgroup CapsuleVsSphere
     * @brief %API for %Capsule vs. %Sphere interactions.
     * 
     * This module contains all helpers for %Capsule vs. %Sphere interactions.
     * @{
     */

    /** 
     * @brief Computes the distance between a capsule and a sphere.
     * 
     * @param c_a The capsule
     * @param s_b The sphere
     * @return The distance between both primitives.
     */
    static double compute_D(const Capsule& c_a, const Sphere& s_b);

    /**
     * @brief Computes the *full* derivative of the distance between a capsule and a sphere.
     * 
     * The derivative is size 9 by 9,
     * meaning it is taken with respect to the parameters of both primitives (all positions, a total of 9 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for both primitives.
     * @param c_a The capsule
     * @param s_b The sphere
     */
    static void compute_dDdP(Vector9d& dDdP, const Capsule& c_a,
                             const Sphere& s_b);

    /**
     * @brief Computes the *partial* derivative of the distance between a capsule and a sphere.
     * 
     * The derivative is of size 6,
     * meaning it is taken with respect to the parameters of the capsule (start- and end-position, a total of 6 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for the capsule.
     * @param c_a The capsule
     * @param s_b The sphere
     * 
     * @attention This is not the *full* derivative, only the segment which belongs to the capsule.
     */
    static void compute_dDdP(Vector6d& dDdP, const Capsule& c_a,
                             const Sphere& s_b);

    /**
     * @brief Computes the *partial* derivative of the distance between a capsule and a sphere.
     * 
     * The derivative is of size 3,
     * meaning it is taken with respect to the parameters of the sphere (position, a total of 3 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for the sphere.
     * @param c_a The capsule
     * @param s_b The sphere
     * 
     * @attention This is not the *full* derivative, only the segment which belongs to the sphere.
     */
    static void compute_dDdP(Vector3d& dDdP, const Capsule& c_a,
                             const Sphere& s_b);

    /**
     * @brief Computes the *full* second derivative of the distance between a capsule and a sphere.
     * 
     * The derivative is of size 9 by 9,
     * meaning it is taken with respect to the parameters of both primitives (all positions, a total of 9 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP2} \f$ where D is the distance and P is the vector of parameters for both primitives.
     * @param c_a The capsule
     * @param s_b The sphere
     */
    static void compute_d2DdP2(Matrix9d& d2DdP2, const Capsule& c_a,
                               const Sphere& s_b);

    /**
     * @brief Computes the *partial* second derivative of the distance between a capsule and a sphere.
     * 
     * The derivative is of size 6 by 6,
     * meaning it is taken with respect to the parameters of the capsule (start- and end-position, a total of 6 degrees of freedom).
     * 
     * @param[out] d2DdP2 The first derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for the capsule.
     * @param c_a The capsule
     * @param s_b The sphere
     * 
     * @attention This is not the *full* second derivative, only the block which belongs to the capsule.
     */
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Capsule& c_a,
                               const Sphere& s_b);

    /**
     * @brief Computes the *partial* second derivative of the distance between a capsule and a sphere.
     * 
     * The derivative is of size 3 by 3,
     * meaning it is taken with respect to the parameters of the sphere (position, a total of 3 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for the sphere.
     * @param c_a The capsule
     * @param s_b The sphere
     * 
     * @attention This is not the *full* derivative, only the block which belongs to the sphere.
     */
    static void compute_d2DdP2(Matrix3d& d2DdP2, const Capsule& c_a,
                               const Sphere& s_b);

    /**
     * @}
     */ // end of group CapsuleVsSphere

    /**
     * @defgroup CapsuleVsCapsule
     * @brief %API for %Capsule vs. %Capsule interactions.
     * 
     * This module contains all helpers for %Capsule vs. %Capsule interactions.
     * @{
     */

    /** 
     * @brief Computes the distance between two capsules.
     * 
     * @param c_a The first capsule
     * @param c_b The second capsule
     * @return The distance between both primitives.
     */
    static double compute_D(const Capsule& c_a, const Capsule& c_b);

    /**
     * @brief Computes the *full* derivative of the distance between both capsules.
     * 
     * The derivative is of size 12,
     * meaning it is taken with respect to the parameters of both capsule (both start- and end-positions, a total of 12 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for both capsules.
     * @param c_a The first capsule
     * @param c_b The second capsule
     */
    static void compute_dDdP(Vector12d& dDdP, const Capsule& c_a,
                             const Capsule& c_b);

    /**
     * @brief Computes the *partial* derivative of the distance between both capsules.
     * 
     * The derivative is of size 6,
     * meaning it is taken with respect to the parameters of the first capsule (namely both start- and end-positions, a total of 6 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for the first capsule.
     * @param c_a The first capsule
     * @param c_b The second capsule
     * 
     * @attention This is not the *full* derivative, only the segment which belong to the first capsule.
     */
    static void compute_dDdP(Vector6d& dDdP, const Capsule& c_a,
                             const Capsule& c_b);

    /**
     * @brief Computes the *full* second derivative of the distance between both capsules.
     * 
     * The derivative is of size 12 by 12,
     * meaning it is taken with respect to the parameters of both capsules (both positions, a total of 6 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for both capsules.
     * @param c_a The first capsule
     * @param c_b The second capsule
     */
    static void compute_d2DdP2(Matrix12d& d2DdP2, const Capsule& c_a,
                               const Capsule& c_b);

    /**
     * @brief Computes the *partial* second derivative of the distance between both capsules.
     * 
     * The derivative is of size 6 by 6,
     * meaning it is taken with respect to the parameters of the first capsule (namely both start- and end-positions, a total of 6 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for the first capsule.
     * @param c_a The first capsule
     * @param c_b The second capsule
     * 
     * @attention This is not the *full* derivative, only the block belonging to the first capsule.
     */
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Capsule& c_a,
                               const Capsule& c_b);

    /**
     * @}
     */ // end of group CapsuleVsCapsule

    /**
     * @defgroup PlaneVsPlane
     * @brief %API for %Plane vs. %Plane interactions.
     * 
     * This module contains all helpers for %Plane vs. %Plane interactions.
     * @{
     */

    /** 
     * @brief Computes the distance between two planes.
     * 
     * @param p_a The first plane
     * @param p_b The second plane
     * @return The distance between both primitives.
     */
    static double compute_D(const Plane& p_a, const Plane& p_b);

    /**
     * @brief Computes the *full* derivative of the distance between two planes.
     * 
     * The derivative is of size 12,
     * meaning it is taken with respect to the parameters of both planes (both positions and normals, a total of 12 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for both planes.
     * @param p_a The first plane
     * @param p_b The second plane
     */
    static void compute_dDdP(Vector12d& dDdP, const Plane& p_a,
                             const Plane& p_b);

    /**
     * @brief Computes the *full* second derivative of the distance between two planes.
     * 
     * The derivative is of size 12 by 12,
     * meaning it is taken with respect to the parameters of both planes (both positions and normals, a total of 12 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for both planes.
     * @param p_a The first plane
     * @param p_b The second plane
     */
    static void compute_d2DdP2(Matrix12d& d2DdP2, const Plane& p_a,
                               const Plane& p_b);

    /**
     * @brief Computes the *partial* derivative of the distance between two planes.
     * 
     * The derivative is of size 6,
     * meaning it is taken with respect to the parameters of the first plane (position and normal, a total of 6 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for the first planes.
     * @param p_a The first plane
     * @param p_b The second plane
     * 
     * @attention This is not the *full* derivative, only the segment which belong to the first plane.
     */
    static void compute_dDdP(Vector6d& dDdP, const Plane& p_a,
                             const Plane& p_b);

    /**
     * @brief Computes the *partial* second derivative of the distance between two planes.
     * 
     * The derivative is of size 6 by 6,
     * meaning it is taken with respect to the parameters of the first plane (position and normal, a total of 6 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for the first planes.
     * @param p_a The first plane
     * @param p_b The second plane
     * 
     * @attention This is not the *full* derivative, only the block which belong to the first plane.
     */
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Plane& p_a,
                               const Plane& p_b);

    /**
     * @}
     */ // end of group PlaneVsPlane

    /**
     * @defgroup PlaneVsCapsule
     * @brief %API for %Plane vs. %Capsule interactions.
     * 
     * This module contains all helpers for %Plane vs. %Capsule interactions.
     * @{
     */

    /** 
     * @brief Computes the distance between a plane and a capsule.
     * 
     * @param p_a The plane
     * @param c_b The capsule
     * @return The distance between both primitives.
     */
    static double compute_D(const Plane& p_a, const Capsule& c_b);

    /**
     * @brief Computes the *full* derivative of the distance between a plane and a capsule.
     * 
     * The derivative is of size 12,
     * meaning it is taken with respect to the parameters of the plane and the capsule (plane position and normal, capsule start- and end-position, a total of 12 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for the plane and the capsule.
     * @param p_a The plane
     * @param c_b The capsule
     */
    static void compute_dDdP(Vector12d& dDdP, const Plane& p_a,
                             const Capsule& c_b);

    /**
     * @brief Computes the *full* second derivative of the distance between a plane and a capsule.
     * 
     * The derivative is of size 12 by 12,
     * meaning it is taken with respect to the parameters of the plane and the capsule (plane position and normal, capsule start- and end-position, a total of 12 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for the plane and the capsule.
     * @param p_a The plane
     * @param c_b The capsule
     */
    static void compute_d2DdP2(Matrix12d& d2DdP2, const Plane& p_a,
                               const Capsule& c_b);

    /**
     * @brief Computes the *partial* derivative of the distance between a plane and a capsule.
     * 
     * The derivative is of size 6,
     * meaning it is taken with respect to the parameters of the plane (position and normal, a total of 6 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for the planes.
     * @param p_a The plane
     * @param c_b The capsule
     * 
     * @attention This is not the *full* derivative, only the segment which belong to the first plane.
     */
    static void compute_dDdP(Vector6d& dDdP, const Plane& p_a,
                             const Capsule& c_b);

    /**
     * @brief Computes the *partial* second derivative of the distance between a plane and a capsule.
     * 
     * The derivative is of size 6 by 6,
     * meaning it is taken with respect to the parameters of the plane (position and normal, a total of 6 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for the planes.
     * @param p_a The plane
     * @param c_b The capsule
     * 
     * @attention This is not the *full* derivative, only the block which belong to the first plane.
     */
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Plane& p_a,
                               const Capsule& c_b);

    /**
     * @}
     */ // end of group PlaneVsCapsule

    /**
     * @defgroup CapsuleVsPlane
     * @brief %API for %Capsule vs. %Plane interactions.
     * 
     * This module contains all helpers for %Capsule vs. %Plane interactions.
     * @{
     */

    /** 
     * @brief Computes the distance between a capsule and a plane.
     * 
     * @param c_a The capsule
     * @param p_b The plane
     * @return The distance between both primitives.
     */
    static double compute_D(const Capsule& c_a, const Plane& p_b);

    /**
     * @brief Computes the *full* derivative of the distance between a capsule and a plane.
     * 
     * The derivative is of size 12,
     * meaning it is taken with respect to the parameters of the capsule and the plane (capsule start- and end-position, plane position and normal, a total of 12 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for the capsule and the plane.
     * @param c_a The capsule
     * @param p_b The plane
     */
    static void compute_dDdP(Vector12d& dDdP, const Capsule& c_a,
                             const Plane& p_b);

    /**
     * @brief Computes the *full* second derivative of the distance between a capsule and a plane.
     * 
     * The derivative is of size 12 by 12,
     * meaning it is taken with respect to the parameters of the capsule and the plane (capsule start- and end-position, plane position and normal, a total of 12 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for capsule and the plane.
     * @param c_a The capsule
     * @param p_b The plane
     */
    static void compute_d2DdP2(Matrix12d& d2DdP2, const Capsule& c_a,
                               const Plane& p_b);

    /**
     * @brief Computes the *partial* derivative of the distance between a capsule and a plane.
     * 
     * The derivative is of size 6,
     * meaning it is taken with respect to the parameters of the capsule (start- and end-position, a total of 6 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for the capsule and the plane.
     * @param c_a The capsule
     * @param p_b The plane
     * 
     * @attention This is not the *full* derivative, only the segment which belong to the capsule.
     */
    static void compute_dDdP(Vector6d& dDdP, const Capsule& c_a,
                             const Plane& p_b);

    /**
     * @brief Computes the *partial* second derivative of the distance between a capsule and a plane.
     * 
     * The derivative is of size 6 by 6,
     * meaning it is taken with respect to the parameters of the capsule (start- and end-position, a total of 6 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for the capsule and the plane.
     * @param c_a The capsule
     * @param p_b The plane
     * 
     * @attention This is not the *full* derivative, only the block which belong to the capsule.
     */
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Capsule& c_a,
                               const Plane& p_b);

    /**
     * @}
     */ // end of group CapsuleVsPlane

    /**
     * @defgroup PlaneVsSphere
     * @brief %API for %Plane vs. %Sphere interactions.
     * 
     * This module contains all helpers for %Plane vs. %Sphere interactions.
     * @{
     */

    /** 
     * @brief Computes the distance between a plane and a sphere.
     * 
     * @param p_a The plane
     * @param s_b The sphere
     * @return The distance between both primitives.
     */
    static double compute_D(const Plane& p_a, const Sphere& s_b);

    /**
     * @brief Computes the *full* derivative of the distance between a plane and a sphere.
     * 
     * The derivative is of size 9,
     * meaning it is taken with respect to the parameters of the plane and the position (plane position and normal, sphere position, a total of 9 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for the plane and the sphere.
     * @param p_a The plane
     * @param s_b The sphere
     */
    static void compute_dDdP(Vector9d& dDdP, const Plane& p_a,
                             const Sphere& s_b);

    /**
     * @brief Computes the *full* second derivative of the distance between a plane and a sphere.
     * 
     * The derivative is of size 9 by 9,
     * meaning it is taken with respect to the parameters of the plane and the sphere (plane position and normal, sphere position, a total of 9 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for the plane and the sphere.
     * @param p_a The plane
     * @param s_b The sphere
     */
    static void compute_d2DdP2(Matrix9d& d2DdP2, const Plane& p_a,
                               const Sphere& s_b);

    /**
     * @brief Computes the *partial* derivative of the distance between a plane and a sphere.
     * 
     * The derivative is of size 6,
     * meaning it is taken with respect to the parameters of the plane (position and normal, a total of 6 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for the plane and the sphere.
     * @param p_a The plane
     * @param s_b The sphere
     * 
     * @attention This is not the *full* derivative, only the segment which belong to the plane.
     */
    static void compute_dDdP(Vector6d& dDdP, const Plane& p_a,
                             const Sphere& s_b);

    /**
     * @brief Computes the *partial* second derivative of the distance between a plane and a sphere.
     * 
     * The derivative is of size 6 by 6,
     * meaning it is taken with respect to the parameters of the plane (position and normal, a total of 6 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for the plane and the sphere.
     * @param p_a The plane
     * @param s_b The sphere
     * 
     * @attention This is not the *full* derivative, only the block which belong to the plane.
     */
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Plane& p_a,
                               const Sphere& s_b);

    /**
     * @brief Computes the *partial* derivative of the distance between a plane and a sphere.
     * 
     * The derivative is of size 3,
     * meaning it is taken with respect to the parameters of the sphere (position, a total of 3 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for the plane and the sphere.
     * @param p_a The plane
     * @param s_b The sphere
     * 
     * @attention This is not the *full* derivative, only the segment which belong to the sphere.
     */
    static void compute_dDdP(Vector3d& dDdP, const Plane& p_a,
                             const Sphere& s_b);

    /**
     * @brief Computes the *partial* second derivative of the distance between a plane and a sphere.
     * 
     * The derivative is of size 3 by 3,
     * meaning it is taken with respect to the parameters of the sphere (position, a total of 3 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for the plane and the sphere.
     * @param p_a The plane
     * @param s_b The sphere
     * 
     * @attention This is not the *full* derivative, only the block which belong to the sphere.
     */
    static void compute_d2DdP2(Matrix3d& d2DdP2, const Plane& p_a,
                               const Sphere& s_b);

    /**
     * @}
     */ // end of group PlaneVsSphere

    /**
     * @defgroup SphereVsPlane
     * @brief %API for %Sphere vs. %Plane interactions.
     * 
     * This module contains all helpers for %Sphere vs. %Plane interactions.
     * @{
     */

    /** 
     * @brief Computes the distance between a sphere and a plane.
     * 
     * @param s_a The sphere
     * @param p_b The plane
     * @return The distance between both primitives.
     */
    static double compute_D(const Sphere& s_a, const Plane& p_b);

    /**
     * @brief Computes the *full* derivative of the distance between a sphere and a plane.
     * 
     * The derivative is of size 9,
     * meaning it is taken with respect to the parameters of the sphere and the plane (sphere position, plane position and normal, a total of 9 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for the sphere and the plane.
     * @param s_a The sphere
     * @param p_b The plane
     */
    static void compute_dDdP(Vector9d& dDdP, const Sphere& s_a,
                             const Plane& p_b);

    /**
     * @brief Computes the *full* second derivative of the distance between a sphere and a plane.
     * 
     * The derivative is of size 9 by 9,
     * meaning it is taken with respect to the parameters of the sphere and the plane (sphere position, plane position and normal, a total of 9 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for sphere and the plane.
     * @param s_a The sphere
     * @param p_b The plane
     */
    static void compute_d2DdP2(Matrix9d& d2DdP2, const Sphere& s_a,
                               const Plane& p_b);

    /**
     * @brief Computes the *partial* derivative of the distance between a sphere and a plane.
     * 
     * The derivative is of size 6,
     * meaning it is taken with respect to the parameters of the plane (position and normal, a total of 6 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for the sphere and the plane.
     * @param s_a The sphere
     * @param p_b The plane
     * 
     * @attention This is not the *full* derivative, only the segment which belong to the plane.
     */
    static void compute_dDdP(Vector6d& dDdP, const Sphere& s_a,
                             const Plane& p_b);

    /**
     * @brief Computes the *partial* second derivative of the distance between a sphere and a plane.
     * 
     * The derivative is of size 6 by 6,
     * meaning it is taken with respect to the parameters of the plane (position and normal, a total of 6 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for the sphere and the plane.
     * @param s_a The sphere
     * @param p_b The plane
     * 
     * @attention This is not the *full* derivative, only the block which belong to the plane.
     */
    static void compute_d2DdP2(Matrix6d& d2DdP2, const Sphere& s_a,
                               const Plane& p_b);

    /**
     * @brief Computes the *partial* derivative of the distance between a sphere and a plane.
     * 
     * The derivative is of size 3,
     * meaning it is taken with respect to the parameters of the sphere (position, a total of 3 degrees of freedom).
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for the sphere and the plane.
     * @param s_a The sphere
     * @param p_b The plane
     * 
     * @attention This is not the *full* derivative, only the segment which belong to the sphere.
     */
    static void compute_dDdP(Vector3d& dDdP, const Sphere& s_a,
                             const Plane& p_b);

    /**
     * @brief Computes the *partial* second derivative of the distance between a sphere and a plane.
     * 
     * The derivative is of size 3 by 3,
     * meaning it is taken with respect to the parameters of the sphere (position, a total of 3 degrees of freedom).
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for the sphere and the plane.
     * @param s_a The sphere
     * @param p_b The plane
     * 
     * @attention This is not the *full* derivative, only the block which belong to the sphere.
     */
    static void compute_d2DdP2(Matrix3d& d2DdP2, const Sphere& s_a,
                               const Plane& p_b);

    /**
     * @}
     */ // end of group SphereVsPlane

    /**
     * @brief Compute the distance for a given pair.
     * 
     * @param[in] p_a The first primitive
     * @param[in] p_b The second primitive
     * @return The distance between both primitives.
     */
    static double compute_D(const primitive_t& p_a, const primitive_t& p_b);

    /**
     * @brief Computes the *full* derivative of the distance for a given pair.
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for both primitives.
     * @param[in] p_a The first primitive
     * @param[in] p_b The second primitive
     */
    static void compute_dDdP(VectorXd& dDdP, const primitive_t& p_a,
                             const primitive_t& p_b);

    /**
     * @brief Computes the *full* second derivative of the distance for a given pair.
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for both primitives.
     * @param[in] p_a The first primitive
     * @param[in] p_b The second primitive
     */
    static void compute_d2DdP2(MatrixXd& d2DdP2, const primitive_t& p_a,
                               const primitive_t& p_b);

    /**
     * @brief Compute the distance for a given pair.
     *  
     * @param[in] pair A pair of indices.
     * @param[in] primitives A vector of primitives, where the indices from the pair should apply to.
     * @return The distance between both primitives.
     */
    static double compute_D(const pair_t& pair,
                            const std::vector<primitive_t>& primitives);

    /**
     * @brief Computes the *full* derivative of the distance for a given pair.
     *
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for both primitives.
     * @param[in] pair A pair of indices.
     * @param[in] primitives A vector of primitives, where the indices from the pair should apply to.
     */
    static void compute_dDdP(VectorXd& dDdP, const pair_t& pair,
                             const std::vector<primitive_t>& primitives);

    /**
     * @brief Computes the *full* second derivative of the distance for a given pair.
     *
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for both primitives.
     * @param[in] pair A pair of indices.
     * @param[in] primitives A vector of primitives, where the indices from the pair should apply to.
     */
    static void compute_d2DdP2(MatrixXd& d2DdP2, const pair_t& pair,
                               const std::vector<primitive_t>& primitives);
};

}  // namespace DCA
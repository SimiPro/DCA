#ifndef __DCA_API_H__
#define __DCA_API_H__

#include "CapsuleVsCapsule.h"
#include "SphereVsCapsule.h"
#include "SphereVsSphere.h"
#include "pair.h"
#include "primitives.h"
#include "utils.h"

namespace DCA {

#define READ_SPHERE_SPHERE_DATA(s_a, s_b) \
    Vector2d props;                       \
    props(0) = s_a.radius;                \
    props(1) = s_b.radius;                \
    Vector6d P;                           \
    P << s_a.position, s_b.position;

#define READ_SPHERE_CAPSULE_DATA(s_a, c_b) \
    Vector2d props;                        \
    props(0) = s_a.radius;                 \
    props(1) = c_b.radius;                 \
    Vector9d P;                            \
    P << s_a.position, c_b.startPosition, c_b.endPosition;

#define READ_CAPSULE_CAPSULE_DATA(c_a, c_b) \
    Vector2d props;                         \
    props(0) = c_a.radius;                  \
    props(1) = c_b.radius;                  \
    Vector12d P;                            \
    P << c_a.startPosition, c_a.endPosition, c_b.startPosition, c_b.endPosition;

class API {
public:
    /**
     * @defgroup SphereVsSphere
     * @brief API for %Sphere vs. %Sphere interactions.
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
    static double compute_D(const Sphere& s_a, const Sphere& s_b) {
        READ_SPHERE_SPHERE_DATA(s_a, s_b);
        return SphereVsSphere::compute_D(P, props);
    }

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
                             const Sphere& s_b) {
        READ_SPHERE_SPHERE_DATA(s_a, s_b);
        SphereVsSphere::compute_dDdP(dDdP, P, props);
    }

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
                             const Sphere& s_b) {
        READ_SPHERE_SPHERE_DATA(s_a, s_b);
        SphereVsSphere::compute_dDdP(dDdP, P, props);
    }

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
                               const Sphere& s_b) {
        READ_SPHERE_SPHERE_DATA(s_a, s_b);
        SphereVsSphere::compute_d2DdP2(d2DdP2, P, props);
    }

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
                               const Sphere& s_b) {
        READ_SPHERE_SPHERE_DATA(s_a, s_b);
        SphereVsSphere::compute_d2DdP2(d2DdP2, P, props);
    }

    /**
     * @}
     */ // end of group SphereVsSphere

    /**
     * @defgroup SphereVsCapsule
     * @brief API for %Sphere vs. %Capsule interactions.
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
    static double compute_D(const Sphere& s_a, const Capsule& c_b) {
        READ_SPHERE_CAPSULE_DATA(s_a, c_b);
        return SphereVsCapsule::compute_D(P, props);
    }

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
                             const Capsule& c_b) {
        READ_SPHERE_CAPSULE_DATA(s_a, c_b);
        SphereVsCapsule::compute_dDdP(dDdP, P, props);
    }

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
                             const Capsule& c_b) {
        READ_SPHERE_CAPSULE_DATA(s_a, c_b);
        SphereVsCapsule::compute_dDdP(dDdP, P, props);
    }

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
                             const Capsule& c_b) {
        READ_SPHERE_CAPSULE_DATA(s_a, c_b);
        SphereVsCapsule::compute_dDdP(dDdP, P, props);
    }

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
                               const Capsule& c_b) {
        READ_SPHERE_CAPSULE_DATA(s_a, c_b);
        SphereVsCapsule::compute_d2DdP2(d2DdP2, P, props);
    }

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
                               const Capsule& c_b) {
        READ_SPHERE_CAPSULE_DATA(s_a, c_b);
        SphereVsCapsule::compute_d2DdP2(d2DdP2, P, props);
    }

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
                               const Capsule& c_b) {
        READ_SPHERE_CAPSULE_DATA(s_a, c_b);
        SphereVsCapsule::compute_d2DdP2(d2DdP2, P, props);
    }

    /**
     * @}
     */ // end of group SphereVsCapsule


    /**
     * @defgroup CapsuleVsSphere
     * @brief API for %Capsule vs. %Sphere interactions.
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
    static double compute_D(const Capsule& c_a, const Sphere& s_b) {
        READ_SPHERE_CAPSULE_DATA(s_b, c_a);
        return SphereVsCapsule::compute_D(P, props);
    }

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
                             const Sphere& s_b) {
        READ_SPHERE_CAPSULE_DATA(s_b, c_a);
        Vector9d dDdP_flip;
        SphereVsCapsule::compute_dDdP(dDdP_flip, P, props);
        dDdP.head(6) = dDdP_flip.tail(6);
        dDdP.tail(3) = dDdP_flip.head(3);
    }

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
                             const Sphere& s_b) {
        READ_SPHERE_CAPSULE_DATA(s_b, c_a);
        SphereVsCapsule::compute_dDdP(dDdP, P, props);
    }

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
                             const Sphere& s_b) {
        READ_SPHERE_CAPSULE_DATA(s_b, c_a);
        SphereVsCapsule::compute_dDdP(dDdP, P, props);
    }

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
                               const Sphere& s_b) {
        READ_SPHERE_CAPSULE_DATA(s_b, c_a);
        Matrix9d d2DdP2_flip;
        SphereVsCapsule::compute_d2DdP2(d2DdP2_flip, P, props);
        d2DdP2.block(0, 0, 6, 6) = d2DdP2_flip.block(3, 3, 6, 6);
        d2DdP2.block(6, 6, 3, 3) = d2DdP2_flip.block(0, 0, 3, 3);

        d2DdP2.block(0, 6, 6, 3) = d2DdP2_flip.block(3, 0, 6, 3);
        d2DdP2.block(6, 0, 3, 6) = d2DdP2_flip.block(0, 3, 3, 6);
    }

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
                               const Sphere& s_b) {
        READ_SPHERE_CAPSULE_DATA(s_b, c_a);
        SphereVsCapsule::compute_d2DdP2(d2DdP2, P, props);
    }

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
                               const Sphere& s_b) {
        READ_SPHERE_CAPSULE_DATA(s_b, c_a);
        SphereVsCapsule::compute_d2DdP2(d2DdP2, P, props);
    }

    /**
     * @}
     */ // end of group CapsuleVsSphere

    /**
     * @defgroup CapsuleVsCapsule
     * @brief API for %Capsule vs. %Capsule interactions.
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
    static double compute_D(const Capsule& c_a, const Capsule& c_b) {
        READ_CAPSULE_CAPSULE_DATA(c_a, c_b);
        return CapsuleVsCapsule::compute_D(P, props);
    }

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
                             const Capsule& c_b) {
        READ_CAPSULE_CAPSULE_DATA(c_a, c_b);
        CapsuleVsCapsule::compute_dDdP(dDdP, P, props);
    }

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
                             const Capsule& c_b) {
        READ_CAPSULE_CAPSULE_DATA(c_a, c_b);
        CapsuleVsCapsule::compute_dDdP(dDdP, P, props);
    }

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
                               const Capsule& c_b) {
        READ_CAPSULE_CAPSULE_DATA(c_a, c_b);
        CapsuleVsCapsule::compute_d2DdP2(d2DdP2, P, props);
    }

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
                               const Capsule& c_b) {
        READ_CAPSULE_CAPSULE_DATA(c_a, c_b);
        CapsuleVsCapsule::compute_d2DdP2(d2DdP2, P, props);
    }

    /**
     * @}
     */ // end of group CapsuleVsCapsule

    /**
     * @brief Compute the distance for a given pair.
     * 
     * @param[in] p_a The first primitive
     * @param[in] p_b The second primitive
     * @return The distance between both primitives.
     */
    static double compute_D(const primitive_t& p_a, const primitive_t& p_b) {
        return std::visit(
            overloaded{[&](const Sphere& s_a, const Sphere& s_b) {
                           return compute_D(s_a, s_b);
                       },
                       [&](const Sphere& s_a, const Capsule& c_b) {
                           return compute_D(s_a, c_b);
                       },
                       [&](const Capsule& c_a, const Sphere& s_b) {
                           return compute_D(c_a, s_b);
                       },
                       [&](const Capsule& c_a, const Capsule& c_b) {
                           return compute_D(c_a, c_b);
                       }},
            p_a, p_b);
    }

    /**
     * @brief Computes the *full* derivative of the distance for a given pair.
     * 
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for both primitives.
     * @param[in] p_a The first primitive
     * @param[in] p_b The second primitive
     */
    static void compute_dDdP(VectorXd& dDdP, const primitive_t& p_a,
                             const primitive_t& p_b) {
        std::visit(overloaded{[&](const Sphere& s_a, const Sphere& s_b) {
                                  Vector6d dDdP_full;
                                  compute_dDdP(dDdP_full, s_a, s_b);
                                  dDdP = dDdP_full;
                              },
                              [&](const Sphere& s_a, const Capsule& c_b) {
                                  Vector9d dDdP_full;
                                  compute_dDdP(dDdP_full, s_a, c_b);
                                  dDdP = dDdP_full;
                              },
                              [&](const Capsule& c_a, const Sphere& s_b) {
                                  Vector9d dDdP_full;
                                  compute_dDdP(dDdP_full, c_a, s_b);
                                  dDdP = dDdP_full;
                              },
                              [&](const Capsule& c_a, const Capsule& c_b) {
                                  Vector12d dDdP_full;
                                  compute_dDdP(dDdP_full, c_a, c_b);
                                  dDdP = dDdP_full;
                              }},
                   p_a, p_b);
    }

    /**
     * @brief Computes the *full* second derivative of the distance for a given pair.
     * 
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for both primitives.
     * @param[in] p_a The first primitive
     * @param[in] p_b The second primitive
     */
    static void compute_d2DdP2(MatrixXd& d2DdP2, const primitive_t& p_a,
                               const primitive_t& p_b) {
        std::visit(overloaded{[&](const Sphere& s_a, const Sphere& s_b) {
                                  Matrix6d d2DdP2_full;
                                  compute_d2DdP2(d2DdP2_full, s_a, s_b);
                                  d2DdP2 = d2DdP2_full;
                              },
                              [&](const Sphere& s_a, const Capsule& c_b) {
                                  Matrix9d d2DdP2_full;
                                  compute_d2DdP2(d2DdP2_full, s_a, c_b);
                                  d2DdP2 = d2DdP2_full;
                              },
                              [&](const Capsule& c_a, const Sphere& s_b) {
                                  Matrix9d d2DdP2_full;
                                  compute_d2DdP2(d2DdP2_full, c_a, s_b);
                                  d2DdP2 = d2DdP2_full;
                              },
                              [&](const Capsule& c_a, const Capsule& c_b) {
                                  Matrix12d d2DdP2_full;
                                  compute_d2DdP2(d2DdP2_full, c_a, c_b);
                                  d2DdP2 = d2DdP2_full;
                              }},
                   p_a, p_b);
    }

    /**
     * @brief Compute the distance for a given pair.
     *  
     * @param[in] pair A pair of indices.
     * @param[in] primitives A vector of primitives, where the indices from the pair should apply to.
     * @return The distance between both primitives.
     */
    static double compute_D(const pair_t& pair,
                            const std::vector<primitive_t>& primitives) {
        return compute_D(primitives.at(pair.first), primitives.at(pair.second));
    }

    /**
     * @brief Computes the *full* derivative of the distance for a given pair.
     *
     * @param[out] dDdP The first derivative \f$ \frac{dD}{dP} \f$ where D is the distance and P is the vector of parameters for both primitives.
     * @param[in] pair A pair of indices.
     * @param[in] primitives A vector of primitives, where the indices from the pair should apply to.
     */
    static void compute_dDdP(VectorXd& dDdP, const pair_t& pair,
                             const std::vector<primitive_t>& primitives) {
        compute_dDdP(dDdP, primitives.at(pair.first),
                     primitives.at(pair.second));
    }

    /**
     * @brief Computes the *full* second derivative of the distance for a given pair.
     *
     * @param[out] d2DdP2 The second derivative \f$ \frac{d^2D}{dP^2} \f$ where D is the distance and P is the vector of parameters for both primitives.
     * @param[in] pair A pair of indices.
     * @param[in] primitives A vector of primitives, where the indices from the pair should apply to.
     */
    static void compute_d2DdP2(MatrixXd& d2DdP2, const pair_t& pair,
                               const std::vector<primitive_t>& primitives) {
        compute_d2DdP2(d2DdP2, primitives.at(pair.first),
                       primitives.at(pair.second));
    }
};

}  // namespace DCA

#endif
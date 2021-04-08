/**
 * This file holds utils used all over this library.
 * 
 * @author: Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_UTILS_H__
#define __DCA_UTILS_H__

#include <eigen3/Eigen/Core>
#include <variant>
#include <vector>

namespace DCA {

// Forward definitions
class Sphere;
class Capsule;

// Type definitions

// All possible primitves
using primitive_t = std::variant<Sphere, Capsule>;

// Easier access to a 3x3 matrix
using Eigen::Matrix3d;

// Easier access to a vector of length 3
using Eigen::Vector3d;

// Easier access to any dynamic matrix
using Eigen::MatrixXd;

// Easier access to any dynamic vector
using Eigen::VectorXd;

// Easier access to a pair (corresponding of two indices)
using pair_t = std::pair<size_t, size_t>;

// Helper for std::visit
template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;


#define EPSILON 1e-8

/**
 * This helper let's us use the fallback to finite difference
 * (for gradient computation) more easily
 */
#define compute_dDdP_FD(OTHER_T)                                    \
    virtual void compute_dDdP(VectorXd &dDdP, const OTHER_T &other) \
        const override {                                            \
        return CollisionWith<OTHER_T>::compute_dDdP(dDdP, other);   \
    }

/**
 * This helper let's us use the fallback to finite difference
 * (for hessian computation) more easily
 */
#define compute_d2DdP2_FD(OTHER_T)                                      \
    virtual void compute_d2DdP2(MatrixXd &d2DdP2, const OTHER_T &other) \
        const override {                                                \
        return CollisionWith<OTHER_T>::compute_d2DdP2(d2DdP2, other);   \
    }

/**
 * This class is the base class for all primitives.
 * Each class must override the distance computation (compute_D),
 * with respect to ALL other primitives.
 * 
 * The derivatives can be calculated using finite difference,
 * by calling compute_dDdP_FD(OTHER_T) with OTHER_T
 * being the type of the other primitive.
 * 
 * The same applies to the hessian.
 */
template <typename T>
class CollisionWith {
public:
    /**
     * Compute the distance to another object.
     * This function must be implemented by any primitive.
     */
    virtual double compute_D(const T &other) const = 0;

    /**
     * Compute the first derivative (gradient)
     * of the distance function to the other object.
     * If this function is not overridden for the different types,
     * it will fall back to finite differences. 
     */
    virtual void compute_dDdP(VectorXd &dDdP, const T &other) const {
        // finite diff
        dDdP.resize(1);
        dDdP(0) = -1;  // for now
    }

    /**
     * Compute the second derivative (hessian)
     * of the distance function to the other object.
     * If this function is not overridden for the different types,
     * it will fall back to finite differences.
     */
    virtual void compute_d2DdP2(MatrixXd &d2DdP2, const T &other) const {
        // finite diff
        d2DdP2.resize(1, 1);
        d2DdP2(0, 0) = -1;  // for now
    }
};

/**
 * This helper class lets you create a new shape
 * without specifying all other types to collide with.
 */
class CollisionWithAll : public CollisionWith<Sphere>,
                         public CollisionWith<Capsule> {};

}  // namespace DCA

#endif /* __DCA_UTILS_H__ */
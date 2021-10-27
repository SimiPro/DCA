#pragma once

#include <DCA/Pair.h>

namespace DCA {

/**
 * @brief Public %API.
 * 
 * This is the main public %API which should be used.
 */
class API {
public:
    /** 
     * @brief Computes the shortest distance between two primitives.
     * 
     * @param[in] p_a The first primitive
     * @param[in] p_b The second primitive
     * @return The shortest distance between both primitives.
     */
    static double compute_D(const primitive_t& p_a, const primitive_t& p_b);

    /** 
     * @brief Computes the first derivative of the shortest distance.
     * 
     * Computes the first derivative of the shortest distance between two primitives with respect to the state of the primitives.
     * 
     * @param dDdD[out] The first derivative \f$ \frac{dD}{dS} \f$
     * @param p_a[in] The first primitive
     * @param p_b[in] The second primitive
     */
    static void compute_dDdS(Vector12d& dDdS, const primitive_t& p_a, const primitive_t& p_b);

    /** 
     * @brief Computes the second derivative of the shortest distance.
     * 
     * Computes the second derivative of the shortest distance between two primitives with respect to the state of the primitives.
     * 
     * @param dDdD[out] The second derivative \f$ \frac{d^2D}{dS^2} \f$
     * @param p_a[in] The first primitive
     * @param p_b[in] The second primitive
     */
    static void compute_d2DdS2(Matrix12d& d2DdS2, const primitive_t& p_a, const primitive_t& p_b);

    /**
     * @brief Compute the points which yield the shortest distance.
     * 
     * Returns a pair of points, where the distance between those points is the shortest between the primitives.
     * The first returned point belongs to primitive p_a, the second belongs to primitive p_b.
     * 
     * @param[in] p_a The first primitive
     * @param[in] p_b The second primitive
     * @return A pair of points, where the distance between those points is the shortest between the two primitives.
     */
    static std::pair<Vector3d, Vector3d> compute_closest_points(const primitive_t& p_a, const primitive_t& p_b);

    /** 
     * @brief Computes the shortest distance between two primitives given \f$ t \f$.
     * 
     * @param[in] p_a The first primitive
     * @param[in] p_b The second primitive
     * @param[in] t The \f$ t \f$ values for the parameterization.
     * @return The shortest distance between both primitives.
     */
    static double compute_D(const primitive_t& p_a, const primitive_t& p_b, const VectorXd& t);

    /** 
     * @brief Computes the first derivative of the shortest distance given \f$ t \f$.
     * 
     * Computes the first derivative of the shortest distance between two primitives with respect to the state of the primitives.
     * 
     * @param dDdD[out] The first derivative \f$ \frac{dD}{dS} \f$
     * @param p_a[in] The first primitive
     * @param p_b[in] The second primitive
     * @param[in] t The \f$ t \f$ values for the parameterization.
     */
    static void compute_dDdS(Vector12d& dDdS, const primitive_t& p_a, const primitive_t& p_b, const VectorXd& t);

    /** 
     * @brief Computes the second derivative of the shortest distance given \f$ t \f$.
     * 
     * Computes the second derivative of the shortest distance between two primitives with respect to the state of the primitives.
     * 
     * @param dDdD[out] The second derivative \f$ \frac{d^2D}{dS^2} \f$
     * @param p_a[in] The first primitive
     * @param p_b[in] The second primitive
     * @param[in] t The \f$ t \f$ values for the parameterization.
     */
    static void compute_d2DdS2(Matrix12d& d2DdS2, const primitive_t& p_a, const primitive_t& p_b, const VectorXd& t);

    /**
     * @brief Compute the points which yield the shortest distance given \f$ t \f$.
     * 
     * Returns a pair of points, where the distance between those points is the shortest between the primitives.
     * The first returned point belongs to primitive p_a, the second belongs to primitive p_b.
     * 
     * @param[in] p_a The first primitive
     * @param[in] p_b The second primitive
     * @param[in] t The \f$ t \f$ values for the parameterization.
     * @return A pair of points, where the distance between those points is the shortest between the two primitives.
     */
    static std::pair<Vector3d, Vector3d> compute_closest_points(const primitive_t& p_a, const primitive_t& p_b, const VectorXd& t);

    /**
     * @brief Compute the \f$ t \f$ values based on the primitives (and their state).
     * 
     * Computes the \f$ t \f$ values based on the state of the primitives.
     * The values represent the points which will yield the shortest distance between the two primitives.
     * 
     * @param[out] t The \f$ t \f$ values which represent the parameterization of the shortest distance.
     * @param[in] p_a The first primitive
     * @param[in] p_b The second primitive
     */
    static void compute_t(VectorXd& t, const primitive_t& p_a, const primitive_t& p_b);
};

}  // namespace DCA
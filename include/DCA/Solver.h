#pragma once

#include <DCA/NewtonOptimizer.h>
#include <DCA/Objective.h>

namespace DCA {

/**
 * @brief Wrapper around a newton optimizer and an objective.
 * 
 * This class actually finds the shortest distance.
 */
class Solver : public FiniteDifference {
public:
    /**
     * @brief Create a solver.
     * @param[in] primitive_A The first primitive.
     * @param[in] primitive_B The second primitive.
     */
    Solver(primitive_t primitive_A, primitive_t primitive_B);

    /**
     * @brief Default deconstructor.
     */
    ~Solver() = default;

    /**
     * @brief Compute the parameterization.
     * @param[out] t The computed parameterization, stacked.
     * @param[in] s The state of the two primitives, stacked.
     * @param[in] forFD Should be used with "false". True will lower residual, which makes it slower but more accurate.
     */
    void compute_t(VectorXd& t, const VectorXd& s, bool forFD = false) const;

    /**
     * @brief Compute the derivative.
     * 
     * Computes \f$ \frac{dt}{ds} \f$.
     * @param[out] dtds \f$ \frac{dt}{ds} \f$.
     * @param[in] s The state of both primitives, stacked.
     * @param[in] t The parameterization of both primitives, stacked.
     */
    void compute_dtds(MatrixXd& dtds, const VectorXd& s, const VectorXd& t) const;

    /**
     * @brief Compute the distance between two primitives.
     * @param[in] s The state of the two primitives.
     * @param[in] t The parameterization of the two primitives, stacked.
     * @return The shortest distance between both primitives, stacked.
     */
    double compute_D(const VectorXd& s, const VectorXd& t) const;

    /**
     * @brief Compute the derivative of the distance between two primitives w.r.t. s.
     * @param[out] dDdS The derivative \f$ \frac{dD}{ds} \f$.
     * @param[in] s The state of the two primitives, stacked.
     * @param[in] t The parameterization of the two primitives, stacked.
     */
    void compute_dDdS(VectorXd& dDdS, const VectorXd& s, const VectorXd& t) const;

    /**
     * @brief Compute the second derivative of the distance between two primitives w.r.t. s.
     * @param[out] d2DdS2 The derivative \f$ \frac{d^2D}{ds^2} \f$.
     * @param[in] s The state of the two primitives, stacked.
     * @param[in] t The parameterization of the two primitives, stacked.
     */
    void compute_d2DdS2(MatrixXd& d2DdS2, const VectorXd& s, const VectorXd& t) const;

    /**
     * @brief Compute the closest points which yield the shortest distance.
     * @param[in] s The state of the two primitives, stacked.
     * @param[in] t The parameterization of the two primitives, stacked.
     * @return A pair of points, where the distance between those points is the shortest between the two primitives.
     */
    std::pair<Vector3d, Vector3d> compute_closest_points(const VectorXd& s, const VectorXd& t) const;

    /**
     * @brief Check if t and s are in sync.
     * @param[in] s The state of the two primitives, stacked.
     * @param[in] t The parameterization of the two primitives, stacked.
     * @return True, if t and s are in sync. False, otherwise.
     */
    bool are_t_and_s_in_sync(const VectorXd& s, const VectorXd& t) const;

    /**
     * @brief Get the state of both primitives, stacked.
     * @return The state of the two primitives, stacked.
     */
    Vector12d get_s_from_primitives() const;

    /**
     * @brief Get the total size of t.
     * @return The total size of t.
     */
    int get_total_size_t() const;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    // helper to test all derivatives.
    void test_derivatives() const;

    // FD tests
    void test_dtds_WithFD(const VectorXd& s, const VectorXd& t) const;
    void test_dDdS_WithFD(const VectorXd& s, const VectorXd& t) const;
    void test_d2DdS2_WithFD(const VectorXd& s, const VectorXd& t) const;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

private:
    Objective objective;                ///< The objective to optimize.
    mutable NewtonOptimizer optimizer;  ///< The newton optimizer.
};

}  // namespace DCA
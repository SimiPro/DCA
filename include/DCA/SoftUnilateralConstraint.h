#pragma once

namespace DCA {

/**
 * @brief This class is a soft constraint, meaning a constraint is softified.
 * 
 * used to model unilateral constraints of the type x < u using a C2 penalty energy f(x).
 * - u is the upper limit that x needs to be less than
 * - if x < u, then the energy of the constraint, its gradient and hessian are all 0 (i.e. inactive)
 * - epsilon is the value away from the limit (how much smaller should x be compared to u) after which f(x) = 0
 * - stiffness controls the rate at which f(x) increases if x > u
 */
class SoftUnilateralConstraint {
public:
    /**
     * @brief Create a constraint.
     * @param[in] limit The limit of the constraint.
     * @param[in] stiffness The stiffness of the softifying.
     * @param[in] epsilon The limit (how close to the limit should the value be able to be).
     */
    SoftUnilateralConstraint(double limit, double stiffness, double epsilon);

    /**
     * @brief Default deconstructor.
     */
    ~SoftUnilateralConstraint() = default;

    /**
     * @brief Set a new limit.
     * @param[in] limit The new limit.
     */
    void setLimit(double limit);

    /**
     * @brief Set a new epsilon.
     * @param[in] eps The new epsilon.
     */
    void setEpsilon(double eps);

    /**
     * @brief Set a new stiffness.
     * @param[in] s The new stiffness.
     */
    void setStiffness(double s);

    /**
     * @brief Compute the force acting on x.
     * @param[in] x The current position
     * @return The force F acting on x.
     */
    double evaluate(double x) const;

    /**
     * @brief Compute first derivative of the force acting on x.
     * @param[in] x The current position
     * @return The first derivative of F evaluated at x.
     */
    double computeDerivative(double x) const;

    /**
     * @brief Compute the second derivative of the force acting on x.
     * @param[in] x The current position
     * @return The second derivative of F evaluated at x.
     */
    double computeSecondDerivative(double x) const;

private:
    double a1, b1, c1, a2, b2, c2, d2, epsilon, limit;  ///< Internal values.
};

}  // namespace DCA
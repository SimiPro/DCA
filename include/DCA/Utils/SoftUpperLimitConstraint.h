#pragma once

namespace DCA {

/**
 * @brief This class is a soft constraint, meaning a constraint is softified.
 */
class SoftUpperLimitConstraint {
private:
    double a1, b1, c1, a2, b2, c2, d2, epsilon;
    double limit = 0;

public:
    /**
     * @brief Create a constraint.
     * @param[in] l The limit of the constraint.
     * @param[in] stiffness The stiffness of the softifying.
     * @param[in] epsilon The limit (how close to the limit should the value be able to be).
     */
    SoftUpperLimitConstraint(double l, double stiffness, double epsilon);

    /**
     * @brief Deconstructor.
     * 
     * It deconstructs things.
     */
    virtual ~SoftUpperLimitConstraint();

    /**
     * @brief Set a new limit.
     * @param[in] l The new limit.
     */
    void setLimit(double l);
    /**
     * @brief Set a new epsilon.
     * @param[in] eps The new epsilon.
     */
    void setEpsilon(double eps);

    /**
     * @brief Compute the force acting on x.
     * @param[in] x The current position
     * @return The force F acting on x.
     */
    double compute_F(double x) const;

    /**
     * @brief Compute first derivative of the force acting on x.
     * @param[in] x The current position
     * @return The first derivative of F evaluated at x.
     */
    double compute_dFdX(double x) const;

    /**
     * @brief Compute the second derivative of the force acting on x.
     * @param[in] x The current position
     * @return The second derivative of F evaluated at x.
     */
    double compute_d2FdX2(double x) const;
};

}  // namespace DCA

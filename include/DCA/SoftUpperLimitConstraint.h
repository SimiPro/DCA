#ifndef __DCA_SOFTUPPERLIMITCONSTRAINT_H__
#define __DCA_SOFTUPPERLIMITCONSTRAINT_H__

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
    SoftUpperLimitConstraint(double l, double stiffness, double epsilon) {
        this->limit = l;
        this->epsilon = epsilon;
        a1 = stiffness;

        b1 = 0.5 * a1 * epsilon;
        c1 = 1. / 6. * a1 * epsilon * epsilon;
        a2 = 1. / (2. * epsilon) * a1;
        b2 = a1;
        c2 = 0.5 * a1 * epsilon;
        d2 = 1. / 6. * a1 * epsilon * epsilon;
    }

    /**
     * @brief Deconstructor.
     * 
     * It deconstructs things.
     */
    virtual ~SoftUpperLimitConstraint() {}

    /**
     * @brief Set a new limit.
     * @param[in] l The new limit.
     */
    void setLimit(double l) { limit = l; }
    /**
     * @brief Set a new epsilon.
     * @param[in] eps The new epsilon.
     */
    void setEpsilon(double eps) { epsilon = eps; }

    /**
     * @brief Compute the force acting on x.
     * @param[in] x The current position
     * @return The force F acting on x.
     */
    double compute_F(double x) const {
        x = x - limit;
        if (x > 0) return 0.5 * a1 * x * x + b1 * x + c1;
        if (x > -epsilon)
            return 1.0 / 3 * a2 * x * x * x + 0.5 * b2 * x * x + c2 * x + d2;
        return 0;
    }

    /**
     * @brief Compute first derivative of the force acting on x.
     * @param[in] x The current position
     * @return The first derivative of F evaluated at x.
     */
    double compute_dFdX(double x) const {
        x = x - limit;
        if (x > 0) return a1 * x + b1;
        if (x > -epsilon) return a2 * x * x + b2 * x + c2;
        return 0;
    }

    /**
     * @brief Compute the second derivative of the force acting on x.
     * @param[in] x The current position
     * @return The second derivative of F evaluated at x.
     */
    double compute_d2FdX2(double x) const {
        x = x - limit;
        if (x > 0) return a1;
        if (x > -epsilon) return 2 * a2 * x + b2;
        return 0;
    }
};

}  // namespace DCA

#endif /* __DCA_SOFTUPPERLIMITCONSTRAINT_H__ */
#ifndef __DCA__SENSITIVIY_H__
#define __DCA__SENSITIVIY_H__

#include "utils.h"

namespace DCA {

/**
 * This class represents an objective used in Sensitivity Analysis.
 */
class SensitivityObjective {
    /**
     * Computes the distance for a given P and X.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual double compute_D(const VectorXd& P, const VectorXd& X) const = 0;

    /**
     * Computes the first derivative of the distance with respect to X.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_dDdX(VectorXd& dDdX, const VectorXd& P,
                              const VectorXd& X) const = 0;

    /**
     * Computes the second derivative of the distance with respect to X.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_d2DdX2(MatrixXd& d2DdX2, const VectorXd& P,
                                const VectorXd& X) const = 0;

    /**
     * Computes the first derivative of the distance with respect to P.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_dDdP(VectorXd& dDdP, const VectorXd& P,
                              const VectorXd& X) const = 0;
    /**
     * Computes the second derivative of the distance with respect to P.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_d2DdP2(MatrixXd& d2DdP2, const VectorXd& P,
                                const VectorXd& X) const = 0;
    /**
     * Computes the second derivative of the distance with respect to X and P.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_d2DdXdP(MatrixXd& d2DdXdP, const VectorXd& P,
                                 const VectorXd& X) const = 0;
};

}  // namespace DCA
#endif /* __DCA__SENSITIVIY_H__ */
#ifndef __DCA__SENSITIVIY_H__
#define __DCA__SENSITIVIY_H__

#include "utils.h"

namespace DCA {

/**
 * This class represents an objective used in Sensitivity Analysis.
 */
template <int SizeP, int SizeX>
class SensitivityObjective {
    using P_v = Eigen::Matrix<double, SizeP, 1>;
    using X_v = Eigen::Matrix<double, SizeX, 1>;
    using X_m = Eigen::Matrix<double, SizeX, SizeX>;
    using P_m = Eigen::Matrix<double, SizeP, SizeP>;
    using XP_m = Eigen::Matrix<double, SizeX, SizeP>;
    /**
     * Computes the distance for a given P and X.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual double compute_D(const P_v& P, const X_v& X) const = 0;

    /**
     * Computes the first derivative of the distance with respect to X.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_dDdX(X_v& dDdX, const P_v& P, const X_v& X) const = 0;

    /**
     * Computes the second derivative of the distance with respect to X.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_d2DdX2(X_m& d2DdX2, const P_v& P,
                                const X_v& X) const = 0;

    /**
     * Computes the first derivative of the distance with respect to P.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_dDdP(P_v& dDdP, const P_v& P, const X_v& X) const = 0;
    /**
     * Computes the second derivative of the distance with respect to P.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_d2DdP2(P_m& d2DdP2, const P_v& P,
                                const X_v& X) const = 0;
    /**
     * Computes the second derivative of the distance with respect to X and P.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_d2DdXdP(XP_m& d2DdXdP, const P_v& P,
                                 const X_v& X) const = 0;
};

}  // namespace DCA
#endif /* __DCA__SENSITIVIY_H__ */
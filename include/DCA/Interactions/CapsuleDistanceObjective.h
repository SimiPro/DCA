#pragma once

#include <DCA/Utils/FD.h>
#include <DCA/Utils/Newton.h>
#include <DCA/Utils/Sensitivity.h>
#include <DCA/Utils/SoftUpperLimitConstraint.h>

namespace DCA {

/**
 * This objective is used to solve the capsule-distance problems:
 * A capsule is represented as a line between two points.
 * We search now the shortest distance between two capsules, i.e. two lines.
 * We solve it using Newont's Method (hence the inheritance of Objective).
 * Furthermore, we provide all derivatives for Sensitivity Analysis.
 */
class CapsuleDistanceObjective : public NewtonObjective<12, 2>,
                                 public SensitivityObjective<12, 2> {
public:
    /**
     * @brief Computes the objective function value of this objective.
     * 
     * @param[in] P The four points (start and end of the first capsule, then start and end position of the second capsule), stacked.
     * @param[in] X The current point on the two lines.
     * 
     * @return The objective function value
     */
    double compute_O(const Vector12d& P, const Vector2d& X) const override;

    FD_CHECK_dOdX(2, 2, 12, 0, "CapsuleDistanceObjective - dOdX");

    /**
     * @brief Computes the first derivative of the objective value.
     * 
     * The derivative is taken with respect to X, the current point on the two lines.
     * 
     * @param[out] dOdX The derivative \f$\frac{dO}{dX}\f$.
     * @param[in] P The four points (start and end of the first capsule, then start and end position of the second capsule), stacked.
     * @param[in] X The current point on the two lines.
     */
    void compute_dOdX(Vector2d& dOdX, const Vector12d& P,
                      const Vector2d& X) const;

    FD_CHECK_d2OdX2(2, 2, 12, 0, "CapsuleDistanceObjective - d2OdX2");

    /**
     * @brief Computes the second derivative of the objective value.
     * @param[out] d2OdX2 The second derivative \f$\frac{d^2O}{dX^2}\f$.
     * @param[in] P The four points (start and end of the first capsule, then start and end position of the second capsule), stacked.
     * @param[in] X The current point on the two lines.
     */
    void compute_d2OdX2(Matrix2d& d2OdX2, const Vector12d& P,
                        const Vector2d& X) const;

    FD_CHECK_dDdX(2, 2, 12, 0, "CapsuleDistanceObjective - dDdX");

    /**
     * @brief Computes the first derivative of the capsule distance with respect to X.
     * 
     * The derivative is taken with respect to X, the current point on the two lines.
     * 
     * @param[out] dDdX The derivative \f$\frac{dD}{dX}\f$.
     * @param[in] P The four points (start and end of the first capsule, then start and end position of the second capsule), stacked.
     * @param[in] X The current point on the two lines.
     */
    void compute_dDdX(Vector2d& dDdX, const Vector12d& P,
                      const Vector2d& X) const;

    FD_CHECK_d2DdX2(2, 2, 12, 0, "CapsuleDistanceObjective - d2DdX2");

    /**
     * @brief Computes the second derivative of the capsule distance with respect to X.
     * 
     * The derivative is taken with respect to X, the current point on the two lines.
     * 
     * @param[out] d2DdX2 The derivative \f$\frac{d^2D}{dX^2}\f$.
     * @param[in] P The four points (start and end of the first capsule, then start and end position of the second capsule), stacked.
     * @param[in] X The current point on the two lines.
     */
    void compute_d2DdX2(Matrix2d& d2DdX2, const Vector12d& P,
                        const Vector2d& X) const override;

    /**
     * @brief Computes the capsule distance between the two capsules parameterized by P.
     * 
     * @param[in] P The four points (start and end of the first capsule, then start and end position of the second capsule), stacked.
     * @param[in] X The current point on the two lines.
     * @return The distance between the two capsules evaluated at X.
     */
    double compute_D(const Vector12d& P, const Vector2d& X) const override;

    FD_CHECK_dDdP_NON_STATIC(12, 12, 2, 0, "CapsuleDistanceObjective - dDdP");

    /**
     * @brief Computes the first derivative of the capsule distance with respect to P.
     * 
     * @param[out] dDdP The derivative \f$\frac{dD}{dP}\f$.
     * @param[in] P The four points (start and end of the first capsule, then start and end position of the second capsule), stacked.
     * @param[in] X The current point on the two lines.
     */
    void compute_dDdP(Vector12d& dDdP, const Vector12d& P,
                      const Vector2d& X) const override;

    FD_CHECK_d2DdP2_NON_STATIC(12, 12, 2, 0,
                               "CapsuleDistanceObjective - d2DdP2");

    /**
     * @brief Computes the second derivative of the capsule distance with respect to P.
     * 
     * @param[out] d2DdP2 The derivative \f$\frac{d^2D}{dP^2}\f$.
     * @param[in] P The four points (start and end of the first capsule, then start and end position of the second capsule), stacked.
     * @param[in] X The current point on the two lines.
     */
    void compute_d2DdP2(Matrix12d& d2DdP2, const Vector12d& P,
                        const Vector2d& X) const override;

    FD_CHECK_d2DdXdP(12, 12, 2, 0, "CapsuleDistanceObjective - d2DdXdP");

    /**
     * @brief Computes the mixed derivative of the capsule distance with respect to X and P.
     * 
     * @param[out] d2DdXdP The derivative \f$\frac{d^D}{dXdP}\f$.
     * @param[in] P The four points (start and end of the first capsule, then start and end position of the second capsule), stacked.
     * @param[in] X The current point on the two lines.
     */
    void compute_d2DdXdP(Eigen::Matrix<double, 2, 12>& d2DdXdP,
                         const Vector12d& P, const Vector2d& X) const override;

private:
    double regularizerWeight = 0.1;  ///< weight for the regularizer
    double constraintWeight = 10.0;  ///< weight for the constraint

    /// 2: Both are regularized / 1: only one is regularized
    int maxRegularizerIndex = 2;

    /// The constraint with constrains the value
    SoftUpperLimitConstraint sulc = SoftUpperLimitConstraint(0.0, 1.0, 0.001);
};

}  // namespace DCA
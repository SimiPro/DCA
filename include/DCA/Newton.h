#ifndef __DCA_NEWTON_H__
#define __DCA_NEWTON_H__

#include <Eigen/Dense>
#include <cmath>  // for std::isfinite
#include <iostream>

#include "utils.h"

namespace DCA {

/**
 * @brief This class represents an objective, which can be used together with the NewtonOptimizer.
 * @param SizeP The size of the parameter vector.
 * @param SizeX The size of the solution vector.
 */
template <int SizeP, int SizeX>
class NewtonObjective {
public:
    using P_v = Eigen::Matrix<double, SizeP, 1>; ///< Helper: Parameter vector
    using X_v = Eigen::Matrix<double, SizeX, 1>; ///< Helper: Solution vector
    using X_m = Eigen::Matrix<double, SizeX, SizeX>; ///< Helper: Solution Matrix (derivative)

    /**
     * @brief Computes the Objective value of this.
     * @param[in] P Some parameters for the objective.
     * @param[in] X The solving variable.
     * @return The objective value of this function.
     */
    virtual double compute_O(const P_v& P, const X_v& X) const = 0;

    /**
     * @brief Computes the derivative of the objective value with respect to X value of this.
     * @param[out] dOdX The first derivative of the objective function value with respect to X.
     * @param[in] P Some parameters for the objective.
     * @param[in] X The solving variable.
     */

    virtual void compute_dOdX(X_v& dOdX, const P_v& P, const X_v& X) const = 0;

    /**
     * @brief Computes the second derivative of the objective value with respect to X value of this.
     * @param[out] d2OdX2 The second derivative of the objective function value with respect to X.
     * @param[in] P Some parameters for the objective.
     * @param[in] X The solving variable.
     */
    virtual void compute_d2OdX2(X_m& d2OdX2, const P_v& P,
                                const X_v& X) const = 0;

    /**
     * @brief This function is called before the optimization step.
     * 
     * It allows to set things up needed later.
     * 
     * @param[in] P The current P value.
     * @param[in] X The current X value.
     */
    virtual void preOptimizationStep(const P_v& P, const X_v& X) {}

    /**
     * @brief This function is called after the optimization step.
     * 
     * It allows to set things up needed later.
     * 
     * @param[in] P The current P value.
     * @param[in] X The current X value.
     */
    virtual void postOptimizationStep(const P_v& P, const X_v& X) {}

    /**
     * @brief This function is called after each linesearch step.
     * 
     * It allows to set things up needed later.
     * 
     * @param[in] P The current P value.
     * @param[in] X The current X value.
     */
    virtual void postLineSearchStep(const P_v& P, const X_v& X) {}

public:
    double weight;
};

/**
 * This class is a simple newton solver, which can be used together with the NewtonObjective.
 * @param SizeP The size of the parameter vector.
 * @param SizeX The size of the solution vector.
 */
template <int SizeP, int SizeX>
class NewtonOptimizer {
public:
    using P_v = Eigen::Matrix<double, SizeP, 1>; ///< Helper: Parameter vector
    using X_v = Eigen::Matrix<double, SizeX, 1>; ///< Helper: Solution vector
    using X_m = Eigen::Matrix<double, SizeX, SizeX>; ///< Helper: Solution Matrix (derivative)

    /**
     * @brief Create a %NewtonOptimizer
     * 
     * @param[in] solverResidual The residual, when the solver should stop early.
     * @param[in] maxLineSearchIterations How many line search iterations should be performed.
     */
    NewtonOptimizer(double solverResidual = 1e-5,
                    unsigned int maxLineSearchIterations = 15)
        : m_solverResidual(solverResidual),
          m_maxLineSearchIterations(maxLineSearchIterations) {}

    /**
     * @brief Deconstructor.
     * 
     * It deconstructs things.
     */
    ~NewtonOptimizer() {}

    /**
     * @brief Optimize the objective.
     * @param[in,out] objective The objective to use. This is changed, since the pre-, postOptimization and postLineSearch functions are called.
     * @param[in] P The parameters.
     * @param[in,out] x The current (and next) solution.
     * @param[in] maxIterations The number of iterations to perform.
     * @return true if the solver converged, false otherwise.
     */
    bool optimize(NewtonObjective<SizeP, SizeX>& objective, const P_v& P,
                  X_v& x, unsigned int maxIterations = 100) {
        m_x_tmp = x;
        bool betterSolutionFound = false;
        bool converged = false;

        for (unsigned int i = 0; i < maxIterations; i++) {
            objective.preOptimizationStep(P, m_x_tmp);

            computeSearchDirection(P, objective);
            betterSolutionFound = doLineSearch(P, objective);

            objective.postOptimizationStep(P, m_x_tmp);

            if (m_gradient.norm() < m_solverResidual) {
                converged = true;
                break;
            }
        }

        if (betterSolutionFound) x = m_x_tmp;
        return converged;
    }

private:
    /**
     * @brief Compute the search direction by searching for the gradient direction
     * @param[in] P The parameters.
     * @param[in,out] objective The objective to compute the search direction for.
     */
    void computeSearchDirection(const P_v& P,
                                NewtonObjective<SizeP, SizeX>& objective) {
        objective.compute_dOdX(m_gradient, P, m_x_tmp);
        objective.compute_d2OdX2(m_hessian, P, m_x_tmp);

        solveLinearSystem(m_searchDir, m_hessian, m_gradient);

        if (m_useDynamicRegularization)
            applyDynamicRegularization(m_searchDir, m_hessian, m_gradient);
    }

    /**
     * @brief Perform line search on the objective
     * @param[in] P The parameters.
     * @param[in,out] objective The objective to do the linesearch on.
     */
    bool doLineSearch(const P_v& P, NewtonObjective<SizeP, SizeX>& objective) {
        if (m_maxLineSearchIterations < 1) {
            m_x_tmp = m_x_tmp - m_searchDir * m_lineSearchStartValue;
            return true;
        }

        double alpha = m_lineSearchStartValue;
        X_v xc(m_x_tmp);
        double initialObjectiveValue = objective.compute_O(P, xc);

        for (unsigned int j = 0; j < m_maxLineSearchIterations; j++) {
            m_x_tmp = xc - m_searchDir * alpha;
            objective.postLineSearchStep(P, m_x_tmp);
            m_objectiveValue = objective.compute_O(P, m_x_tmp);

            if (!std::isfinite(m_objectiveValue))
                m_objectiveValue = initialObjectiveValue + 1.;

            if (m_objectiveValue > initialObjectiveValue)
                alpha *= 0.5;
            else
                return true;
        }
        return false;
    }

    /**
     * @brief Solve the system A * y = x for y.
     * @param[out] y The solution of the system.
     * @param[in] A The matrix A.
     * @param[in] x The right-hand-side of the system.
     */
    static void solveLinearSystem(X_v& y, const X_m& A, const X_v& x) {
        y = A.colPivHouseholderQr().solve(x);
    }

    /**
     * @brief Apply dynamic regularization on the linear system A * y = x
     * @param[in,out] y The solution of the system
     * @param[in,out] A The matrix A.
     * @param[in] x The right-hand-side of the system.
     */
    static void applyDynamicRegularization(X_v& y, X_m& A, const X_v& x) {
        double dotProduct = y.dot(x);
        if (dotProduct <= 0. && x.squaredNorm() > 0) {
            VectorXd stabRegularizer(A.rows());
            stabRegularizer.setZero();

            double currStabValue = 1e-4;
            for (int i = 0; i < 10; i++) {
                stabRegularizer.setConstant(currStabValue);
                A += stabRegularizer.asDiagonal();
                currStabValue *= 10.;

                solveLinearSystem(y, A, x);

                dotProduct = y.dot(x);
                if (dotProduct > 0.) break;
            }
        }
    }

private:
    ///< residual of the solver
    double m_solverResidual;
    ///< how many line search steps should be done
    unsigned int m_maxLineSearchIterations;
    ///< the starting value of the line search
    double m_lineSearchStartValue = 1.0;
    ///< whether to use dynamic regularization
    bool m_useDynamicRegularization = true;

    ///< the current objective value
    double m_objectiveValue;
    ///< some private members, could also be given via parameters
    X_v m_x_tmp, m_searchDir, m_gradient;
    ///< and the hessian
    X_m m_hessian;
};

}  // namespace DCA

#endif /* __DCA_NEWTON_H__ */
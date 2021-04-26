/**
 * This file holds a simple newton optimizer.
 * 
 * @author Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_NEWTON_H__
#define __DCA_NEWTON_H__

#include <cmath>  // for std::isfinite
#include <Eigen/Dense>
#include <iostream>

#include "utils.h"

namespace DCA {


/**
 * This class represents an objective, which can be used together with the NewtonMinimizer.
 */
class NewtonObjective {
public:
    /**
     * Computes the Objective value of this.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual double compute_O(const VectorXd& P, const VectorXd& X) const = 0;
    /**
     * Computes the derivative of the objective value with respect to X value of this.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_dOdX(VectorXd& dOdX, const VectorXd& P,
                              const VectorXd& X) const = 0;
    /**
     * Computes the second derivative of the objective value with respect to X value of this.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_d2OdX2(MatrixXd& d2OdX2, const VectorXd& P,
                                const VectorXd& X) const = 0;

    virtual void preOptimizationStep(const VectorXd& P, const VectorXd& X) {}
    virtual void postOptimizationStep(const VectorXd& P, const VectorXd& X) {}
    virtual void postLineSearchStep(const VectorXd& P, const VectorXd& X) {}

public:
    double weight;
};

/**
 * @todo: 
 * don't store class members, but
 * rather give by reference/value
 */
class NewtonOptimizer {
public:
    NewtonOptimizer(double solverResidual = 1e-5,
                    unsigned int maxLineSearchIterations = 15)
        : m_solverResidual(solverResidual),
          m_maxLineSearchIterations(maxLineSearchIterations) {}

    ~NewtonOptimizer() {}

    /**
     * Optimize the objective.
     * @return true if the solver converged, false otherwise.
     */
    bool optimize(NewtonObjective& objective, const VectorXd& P, VectorXd& x,
                  unsigned int maxIterations = 100) {
        m_x_tmp = x;
        bool betterSolutionFound = false;
        bool converged = false;

        for (int i = 0; i < maxIterations; i++) {
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
     * Compute the search direction by searching for the gradient direction
     */
    void computeSearchDirection(const VectorXd& P, NewtonObjective& objective) {
        objective.compute_dOdX(m_gradient, P, m_x_tmp);
        objective.compute_d2OdX2(m_hessian, P, m_x_tmp);

        solveLinearSystem(m_searchDir, m_hessian, m_gradient);

        if (m_useDynamicRegularization)
            applyDynamicRegularization(m_searchDir, m_hessian, m_gradient);
    }

    /**
     * Perform line search on the objective
     */
    bool doLineSearch(const VectorXd& P, NewtonObjective& objective) {
        if (m_maxLineSearchIterations < 1) {
            m_x_tmp = m_x_tmp - m_searchDir * m_lineSearchStartValue;
            return true;
        }

        double alpha = m_lineSearchStartValue;
        VectorXd xc(m_x_tmp);
        double initialObjectiveValue = objective.compute_O(P, xc);

        for (int j = 0; j < m_maxLineSearchIterations; j++) {
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
     * Solve the system A * y = x for y.
     */
    static void solveLinearSystem(VectorXd& y, const MatrixXd& A,
                                  const VectorXd& x) {
        y = A.colPivHouseholderQr().solve(x);
    }

    /**
     * Apply dynamic regularization on the linear system A * y = x
     */
    static void applyDynamicRegularization(VectorXd& y, MatrixXd& A,
                                           const VectorXd& x) {
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
    double m_solverResidual; ///< residual of the solver
    unsigned int m_maxLineSearchIterations; ///< how many line search steps should be done
    double m_lineSearchStartValue = 1.0; ///< the starting value of the line search
    bool m_useDynamicRegularization = true; ///< whether to use dynamic regularization

    double m_objectiveValue; ///< the current objective value
    VectorXd m_x_tmp, m_searchDir, m_gradient; ///< some private members, could also be given via parameters
    MatrixXd m_hessian; ///< and the hessian
};

}  // namespace DCA

#endif /* __DCA_NEWTON_H__ */
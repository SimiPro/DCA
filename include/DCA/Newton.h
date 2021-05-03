/**
 * This file holds a simple newton optimizer.
 * 
 * @author Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_NEWTON_H__
#define __DCA_NEWTON_H__

#include <Eigen/Dense>
#include <cmath>  // for std::isfinite
#include <iostream>

#include "utils.h"

namespace DCA {

/**
 * This class represents an objective, which can be used together with the NewtonMinimizer.
 */
template <int SizeP, int SizeX>
class NewtonObjective {
public:
    using P_v = Eigen::Matrix<double, SizeP, 1>;
    using X_v = Eigen::Matrix<double, SizeX, 1>;
    using X_m = Eigen::Matrix<double, SizeX, SizeX>;

    /**
     * Computes the Objective value of this.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual double compute_O(const P_v& P, const X_v& X) const = 0;

    /**
     * Computes the derivative of the objective value with respect to X value of this.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */

    virtual void compute_dOdX(X_v& dOdX, const P_v& P, const X_v& X) const = 0;

    /**
     * Computes the second derivative of the objective value with respect to X value of this.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_d2OdX2(X_m& d2OdX2, const P_v& P,
                                const X_v& X) const = 0;

    virtual void preOptimizationStep(const P_v& P, const X_v& X) {}
    virtual void postOptimizationStep(const P_v& P, const X_v& X) {}
    virtual void postLineSearchStep(const P_v& P, const X_v& X) {}

public:
    double weight;
};

/**
 * @todo: 
 * don't store class members, but
 * rather give by reference/value
 */
template <int SizeP, int SizeX>
class NewtonOptimizer {
public:
    using P_v = Eigen::Matrix<double, SizeP, 1>;
    using X_v = Eigen::Matrix<double, SizeX, 1>;
    using X_m = Eigen::Matrix<double, SizeX, SizeX>;
    NewtonOptimizer(double solverResidual = 1e-5,
                    unsigned int maxLineSearchIterations = 15)
        : m_solverResidual(solverResidual),
          m_maxLineSearchIterations(maxLineSearchIterations) {}

    ~NewtonOptimizer() {}

    /**
     * Optimize the objective.
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
     * Compute the search direction by searching for the gradient direction
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
     * Perform line search on the objective
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
     * Solve the system A * y = x for y.
     */
    static void solveLinearSystem(X_v& y, const X_m& A, const X_v& x) {
        y = A.colPivHouseholderQr().solve(x);
    }

    /**
     * Apply dynamic regularization on the linear system A * y = x
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
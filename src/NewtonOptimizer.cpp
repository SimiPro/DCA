#include <DCA/Logger.h>
#include <DCA/NewtonOptimizer.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace DCA {

bool NewtonOptimizer::optimize(VectorXd& t, const VectorXd& s, const Objective& objective, int maxIterations) {
    //--- Preparations
    t_tmp = t;  //We make a copy here, because we don't know beforehand if the solution we get is better than before
    bool betterSolutionFound = false;
    bool converged = false;

    //--- Apply optimizationLoop
    for (int i = 0; i < maxIterations; i++) {
        //--- Compute search direction
        computeSearchDirection(s, objective);

        //--- Check convergence
        if (gradient.norm() < solverResidual) {
            converged = true;
            //Logger::print(Logger::MAGENTA, "Number of iterations: %d\n", ++i);
            break;
        }

        //--- Apply line search
        if (doLineSearch(s, objective))
            betterSolutionFound = true;
    }

    //--- Print what needs to be printed
    if (printInfo)
        printFinalValues(converged, betterSolutionFound);

    //--- Finish up
    if (betterSolutionFound)  //We only copy it over if we found a better solution
        t = t_tmp;
    return converged;
}

void NewtonOptimizer::computeSearchDirection(const VectorXd& s, const Objective& objective) {
    //--- Compute gradient
    objective.compute_pOpT(gradient, s, t_tmp);

    //--- Compute hessian
    objective.compute_p2OpT2(hessian, s, t_tmp);
    if (checkHessianRank)
        applyMatrixInvertibilityCheck(hessian, name);

    //--- Solve linear system
    solveLinearSystem(searchDir, hessian, gradient, name);
    if (checkLinearSystemSolve)
        applyLinearSystemSolveCheck(searchDir, hessian, gradient, name);

    //--- Apply dynamic regularization
    if (useDynamicRegularization)
        applyDynamicRegularization(searchDir, hessian, gradient, name);
}

bool NewtonOptimizer::doLineSearch(const VectorXd& s, const Objective& objective) {
    //--- Apply no line search
    if (maxLineSearchIterations < 1) {
        Logger::print(Logger::YELLOW, "WARNING: NewtonOptimizer::doLineSearch -> no line search is applied...\n");
        t_tmp = t_tmp - searchDir * lineSearchStartValue;
        return true;
    }

    //--- Initialize things
    double alpha = lineSearchStartValue;
    VectorXd tc(t_tmp);
    double initialObjectiveValue = objective.compute_O(s, tc);

    //--- Line search loop
    for (int j = 0; j < maxLineSearchIterations; j++) {
        //--- Try new solution
        t_tmp = tc - searchDir * alpha;
        objectiveValue = objective.compute_O(s, t_tmp);

        //If not satisfying, try again
        if (objectiveValue > initialObjectiveValue)
            alpha /= 2.0;
        else
            return true;  //Better solution found!
    }

    return false;  //NO better solution found!
}

void NewtonOptimizer::solveLinearSystem(VectorXd& y, const MatrixXd& A, const VectorXd& x, const std::string& description) {
    Eigen::LLT<Eigen::MatrixXd, Eigen::Lower> solver;
    solver.compute(A.triangularView<Eigen::Lower>());
    y = solver.solve(x);
    if (solver.info() != Eigen::Success)
        Logger::print(Logger::RED, "WARNING -> NewtonOptimizer::solveLinearSystem -> %s: y = A*x (dense) solve unsuccessful!\n", description.c_str());
}

bool NewtonOptimizer::applyLinearSystemSolveCheck(VectorXd& y, const MatrixXd& A, const VectorXd& x, const std::string& description) {
    VectorXd y_test = A.inverse() * x;
    double testNorm = (y_test - y).norm();
    bool testPassed = false;
    if (std::isnan(testNorm)) {
        Logger::print(Logger::RED, "WARNING: NewtonOptimizer::evaluateLinearSystemSolveCheck -> %s -> testNorm is NaN: %lf\n", description.c_str(), testNorm);
    } else if (testNorm > 1e-6) {
        Logger::print(Logger::RED, "WARNING: NewtonOptimizer::evaluateLinearSystemSolveCheck -> %s -> system is NOT solved correctly: %lf\n",
                      description.c_str(), testNorm);
    } else {
        Logger::print(Logger::GREEN, "NewtonOptimizer::evaluateLinearSystemSolveCheck -> %s -> system is solved just fine: %lf\n", description.c_str(),
                      testNorm);
        testPassed = true;
    }
    return testPassed;
}

bool NewtonOptimizer::applyMatrixInvertibilityCheck(const MatrixXd& matrix, const std::string& description) {
    Eigen::FullPivLU<MatrixXd> lu(matrix);
    if (!lu.isInvertible()) {
        Logger::print(Logger::RED, "WARNING: NewtonOptimizer::applyMatrixInvertibilityCheck -> %s -> matrix is NOT invertible\n", description.c_str());
        return false;
    }
    Logger::print(Logger::GREEN, "NewtonOptimizer::applyMatrixInvertibilityCheck -> %s -> matrix is READY to be inverted!\n", description.c_str());
    return true;
}

void NewtonOptimizer::applyDynamicRegularization(VectorXd& y, MatrixXd& A, const VectorXd& x, const std::string& description) {
    double dotProduct = y.dot(x);
    if (dotProduct <= 0 && x.squaredNorm() > 0) {
        VectorXd stabRegularizer(A.rows());
        stabRegularizer.setZero();

        double currStabValue = 1e-4;
        for (int i = 0; i < 10; ++i) {
            stabRegularizer.setConstant(currStabValue);
            A += stabRegularizer.asDiagonal();
            currStabValue *= 10;

            solveLinearSystem(y, A, x, description);

            dotProduct = y.dot(x);
            if (dotProduct > 0)
                break;
        }
    }
}

void NewtonOptimizer::printFinalValues(const bool& converged, const bool& betterSolutionFound) {
    if (converged)
        Logger::print(Logger::GREEN, "%s converged ", name.c_str());
    else if (betterSolutionFound)
        Logger::print(Logger::YELLOW, "%s in progress ", name.c_str());
    else
        Logger::print(Logger::RED, "%s is stuck ", name.c_str());
    Logger::print(Logger::DEFAULT, "-> Function value: %10.10lf. Grad norm: %lf. ", objectiveValue, gradient.norm());
    if (searchDir.dot(gradient) <= 0 && gradient.squaredNorm() > 0)
        Logger::print(Logger::RED, "Search dir norm: %lf\n", searchDir.norm());
    else
        Logger::print(Logger::DEFAULT, "Search dir norm: %lf\n", searchDir.norm());
}

}  // namespace DCA
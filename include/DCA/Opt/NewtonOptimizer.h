#pragma once

#include <DCA/Opt/Objective.h>

namespace DCA {
namespace Opt {

class NewtonOptimizer {
public:
    NewtonOptimizer(const std::string& name = "NewtonOptimizer", double solverResidual = 1e-5, int maxLineSearchIterations = 15)
        : name(name), solverResidual(solverResidual), maxLineSearchIterations(maxLineSearchIterations) {}
    ~NewtonOptimizer() {}

    //Returns whether or not solver has converged
    bool optimize(VectorXd& t, const VectorXd& s, const Objective& objective, int maxIterations);

private:
    //Everything to make optization happen
    void computeSearchDirection(const VectorXd& s, const Objective& objective);
    bool doLineSearch(const VectorXd& s, const Objective& objective);
    void applyDynamicRegularization(VectorXd& y, MatrixXd& A, const VectorXd& x, const std::string& description);
    void printFinalValues(const bool& converged, const bool& betterSolutionFound);

    //Helpers to ensure that everything is working correctly
    static void solveLinearSystem(VectorXd& y, const MatrixXd& A, const VectorXd& x, const std::string& description);
    static bool applyLinearSystemSolveCheck(VectorXd& y, const MatrixXd& A, const VectorXd& x, const std::string& description);
    static bool applyMatrixInvertibilityCheck(const MatrixXd& matrix, const std::string& description);

public:
    //Members set by constructor
    std::string name;             //Name of the solver (mainly used for printing)
    double solverResidual;        //Accurancy when solver states that it has converged (compared to the norm of the gradient)
    int maxLineSearchIterations;  //Max number of line search iterations before giving up

    //Other settings
    double lineSearchStartValue = 1.0;     //Start value of line search. If unsuccessful, this value is cut in half until a suitable solution is found
    bool useDynamicRegularization = true;  //If activated, we regularize the hessian in case we don't have a decending search direction
    bool printInfo = false;                //Print infos to console to observe solver progress
    bool checkLinearSystemSolve = false;   //Check if linear system has been solved successfully (i.e. the correct linear solver has been used)
    bool checkHessianRank = false;         //Check if hessian is invertible

private:
    //Helpers for optimization are stored here
    double objectiveValue;
    VectorXd t_tmp, searchDir, gradient;
    MatrixXd hessian;
};

}  // namespace Opt
}  // namespace DCA
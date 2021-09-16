#include <DCA/Opt/FiniteDifference.h>
#include <DCA/Utils/Logger.h>

namespace DCA {
namespace Opt {

void FiniteDifference::setFDCheckIsBeingApplied(bool isBeingApplied) const {
    FDCheckIsBeingApplied = isBeingApplied;
}

void FiniteDifference::estimateVector(VectorXd& vector, const Function_Value& evaluate, const VectorXd& s, const VectorXd& t, D_WRT d_wrt) const {
    VectorXd p_der = (d_wrt == D_S) ? s : t;
    vector.resize(p_der.size());
    vector.setZero();

    auto eval = [&](D_WRT d_wrt) -> double {
        if (d_wrt == D_S)
            return evaluate(p_der, t);
        else if (d_wrt == D_T)
            return evaluate(s, p_der);
        return 0.0;
    };

    double f_P, f_M;
    for (int i = 0; i < p_der.size(); i++) {
        p_der[i] += deltaFD;
        f_P = eval(d_wrt);

        p_der[i] -= 2.0 * deltaFD;
        f_M = eval(d_wrt);

        p_der[i] += deltaFD;
        vector[i] = (f_P - f_M) / (2 * deltaFD);
    }
}

void FiniteDifference::estimateMatrix(MatrixXd& matrix, const Function_Vector& evaluate, const VectorXd& s, const VectorXd& t, D_WRT d_wrt,
                                      int firstDim) const {
    VectorXd p_der = (d_wrt == D_S) ? s : t;
    matrix.resize(firstDim, p_der.size());
    matrix.setZero();

    auto eval = [&](VectorXd& vec, D_WRT d_wrt) -> void {
        if (d_wrt == D_S)
            evaluate(vec, p_der, t);
        else if (d_wrt == D_T)
            evaluate(vec, s, p_der);
    };

    VectorXd f_P(firstDim), f_M(firstDim);
    for (int i = 0; i < p_der.size(); i++) {
        f_P.setZero();
        f_M.setZero();

        p_der[i] += deltaFD;
        eval(f_P, d_wrt);

        p_der[i] -= 2.0 * deltaFD;
        eval(f_M, d_wrt);

        p_der[i] += deltaFD;
        matrix.col(i) = (f_P - f_M) / (2 * deltaFD);
    }
}

void FiniteDifference::testVector(const Function_Value& evaluate, const Function_Vector& analytic, const VectorXd& s, const VectorXd& t, D_WRT d_wrt,
                                  const std::string& name) const {
    VectorXd p_der = (d_wrt == D_S) ? s : t;
    uint pSize = (uint)p_der.size();

    VectorXd fdVector(pSize);
    fdVector.setZero();
    VectorXd anaVector(pSize);
    anaVector.setZero();

    setFDCheckIsBeingApplied(true);
    estimateVector(fdVector, evaluate, s, t, d_wrt);
    analytic(anaVector, s, t);
    setFDCheckIsBeingApplied(false);

    printFDCheck(fdVector, anaVector, name);
}

void FiniteDifference::testMatrix(const Function_Vector& evaluate, const Function_Matrix& analytic, const VectorXd& s, const VectorXd& t, D_WRT d_wrt,
                                  const std::string& name, int firstDim) const {
    VectorXd p_der = (d_wrt == D_S) ? s : t;
    uint pSize = (uint)p_der.size();

    MatrixXd fdMatrix(firstDim, pSize);
    fdMatrix.setZero();
    MatrixXd anaMatrix(firstDim, pSize);
    anaMatrix.setZero();

    setFDCheckIsBeingApplied(true);
    estimateMatrix(fdMatrix, evaluate, s, t, d_wrt, firstDim);
    analytic(anaMatrix, s, t);
    setFDCheckIsBeingApplied(false);

    printFDCheck(fdMatrix, anaMatrix, name);
}

void FiniteDifference::printFDCheck(const VectorXd& estimate, const VectorXd& analytic, const std::string& name) const {
    bool checkSuccessful = true;
    Logger::print(Logger::BLUE,
                  "------------------------------------------------------------"
                  "-----------------------------------------------\n");
    Logger::print(Logger::DEFAULT, "Test ");
    Logger::print(Logger::CYAN, "%s", name.c_str());
    Logger::print(Logger::DEFAULT, " of ");
    Logger::print(Logger::YELLOW, "%s", description.c_str());
    Logger::print(Logger::DEFAULT, " with FD -> norms: Ana: %lf, FD: %lf\n", analytic.norm(), estimate.norm());
    for (int i = 0; i < (int)estimate.size(); i++) {
        double absErr = std::abs(estimate[i] - analytic[i]);
        double relError = 2 * absErr / (eps + std::abs(analytic[i]) + std::abs(estimate[i]));
        if (relError > relTol && absErr > absTol) {
            checkSuccessful = false;
            Logger::print(Logger::RED, "   MISMATCH");
            Logger::print(Logger::DEFAULT,
                          " -> Element: %d, Ana val: %lf, FD val: %lf, Error: "
                          "%lf(%lf%%)\n",
                          i, analytic[i], estimate[i], absErr, relError * 100);
        }
    }
    if (checkSuccessful)
        Logger::print(Logger::GREEN, "   PASSED\n");
    Logger::print(Logger::DEFAULT, "End of check of %s\n", description.c_str());
    Logger::print(Logger::BLUE,
                  "------------------------------------------------------------"
                  "-----------------------------------------------\n");
}

void FiniteDifference::printFDCheck(const MatrixXd& estimate, const MatrixXd& analytic, const std::string& name) const {
    bool checkSuccessful = true;
    Logger::print(Logger::BLUE,
                  "------------------------------------------------------------"
                  "-----------------------------------------------\n");
    Logger::print(Logger::DEFAULT, "Test ");
    Logger::print(Logger::CYAN, "%s", name.c_str());
    Logger::print(Logger::DEFAULT, " of ");
    Logger::print(Logger::YELLOW, "%s", description.c_str());
    if (analytic.size() > 0 && estimate.size() > 0)
        Logger::print(Logger::DEFAULT, " with FD -> norms: Ana: %lf, FD: %lf\n", analytic.norm(), estimate.norm());
    else
        Logger::print(Logger::DEFAULT, " with FD -> norms: Ana: %lf, FD: %lf\n", 0.0, 0.0);
    for (int i = 0; i < (int)estimate.rows(); i++) {
        for (int j = 0; j < (int)estimate.cols(); j++) {
            double absErr = std::abs(estimate.coeff(i, j) - analytic.coeff(i, j));
            double relError = 2 * absErr / (eps + std::abs(estimate.coeff(i, j)) + std::abs(analytic.coeff(i, j)));
            if (relError > relTol && absErr > absTol) {
                checkSuccessful = false;
                Logger::print(Logger::RED, "   MISMATCH");
                Logger::print(Logger::DEFAULT,
                              " -> Element: (%d, %d), Ana val: %lf, FD val: "
                              "%lf. Error: %lf(%lf%%)\n",
                              i, j, analytic.coeff(i, j), estimate.coeff(i, j), absErr, relError * 100);
            }
        }
    }
    if (checkSuccessful)
        Logger::print(Logger::GREEN, "   PASSED\n");
    Logger::print(Logger::DEFAULT, "End of check of %s\n", description.c_str());
    Logger::print(Logger::BLUE,
                  "------------------------------------------------------------"
                  "-----------------------------------------------\n");
}

}  // namespace Opt
}  // namespace DCA
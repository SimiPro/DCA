#pragma once

#include <DCA/Utils/Utils.h>

#include <functional>

namespace DCA {

namespace Opt {

class FiniteDifference {
public:
    FiniteDifference(std::string description) : description(description) {}
    virtual ~FiniteDifference() {}

    virtual void setFDCheckIsBeingApplied(bool isBeingApplied) const;

    //Derivative with respect to ...
    enum D_WRT { D_S, D_T };

    typedef std::function<double(const VectorXd& s, const VectorXd& t)> Function_Value;
    typedef std::function<void(VectorXd& vector, const VectorXd& s, const VectorXd& t)> Function_Vector;
    typedef std::function<void(MatrixXd& matrix, const VectorXd& s, const VectorXd& t)> Function_Matrix;

    void estimateVector(VectorXd& vector, const Function_Value& evaluate, const VectorXd& s, const VectorXd& t, D_WRT d_wrt) const;
    void estimateMatrix(MatrixXd& matrix, const Function_Vector& evaluate, const VectorXd& s, const VectorXd& t, D_WRT d_wrt, int firstDim) const;

    void testVector(const Function_Value& evaluate, const Function_Vector& analytic, const VectorXd& s, const VectorXd& t, D_WRT d_wrt,
                    const std::string& name) const;
    void testMatrix(const Function_Vector& evaluate, const Function_Matrix& analytic, const VectorXd& s, const VectorXd& t, D_WRT d_wrt,
                    const std::string& name, int firstDim) const;

    void printFDCheck(const VectorXd& estimate, const VectorXd& analytic, const std::string& name) const;
    void printFDCheck(const MatrixXd& estimate, const MatrixXd& analytic, const std::string& name) const;

protected:
    mutable bool FDCheckIsBeingApplied = false;
    double deltaFD = 1e-6;
    double relTol = 1e-4;
    double absTol = 1e-6;
    double eps = 1e-10;

public:
    std::string description;
};

}  // namespace Opt
}  // namespace DCA
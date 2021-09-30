#pragma once

namespace DCA {

/*
class used to model unilateral constraints of the type x < u using a C2 penalty energy f(x).
	- u is the upper limit that x needs to be less than
	- if x < u, then the energy of the constraint, its gradient and hessian are all 0 (i.e. inactive)
	- epsilon is the value away from the limit (how much smaller should x be compared to u) after which f(x) = 0
	- stiffness controls the rate at which f(x) increases if x > u
*/
class SoftUnilateralConstraint {
public:
    SoftUnilateralConstraint(double limit, double stiffness, double epsilon);
    ~SoftUnilateralConstraint() {}

    void setLimit(double lim);
    void setEpsilon(double eps);
    void setStiffness(double s);

    //comptue f(x)
    double evaluate(double x) const;

    //compute df/dx
    double computeDerivative(double x) const;

    //compute ddf/dxdx
    double computeSecondDerivative(double x) const;

private:
    double a1, b1, c1, a2, b2, c2, d2, epsilon, limit;
};

}  // namespace DCA
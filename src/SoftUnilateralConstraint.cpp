#include <DCA/Opt/SoftUnilateralConstraint.h>

namespace DCA {
namespace Opt {

SoftUnilateralConstraint::SoftUnilateralConstraint(double limit, double stiffness, double epsilon) {
    this->limit = limit;
    this->epsilon = epsilon;
    a1 = stiffness;

    b1 = 0.5 * a1 * epsilon;
    c1 = 1. / 6. * a1 * epsilon * epsilon;
    a2 = 1. / (2. * epsilon) * a1;
    b2 = a1;
    c2 = 0.5 * a1 * epsilon;
    d2 = 1. / 6. * a1 * epsilon * epsilon;
}

void SoftUnilateralConstraint::setLimit(double lim) {
    limit = lim;
}

void SoftUnilateralConstraint::setEpsilon(double eps) {
    epsilon = eps;

    b1 = 0.5 * a1 * epsilon;
    c1 = 1. / 6. * a1 * epsilon * epsilon;
    a2 = 1. / (2. * epsilon) * a1;
    b2 = a1;
    c2 = 0.5 * a1 * epsilon;
    d2 = 1. / 6. * a1 * epsilon * epsilon;
}

void SoftUnilateralConstraint::setStiffness(double s) {
    a1 = s;

    b1 = 0.5 * a1 * epsilon;
    c1 = 1. / 6. * a1 * epsilon * epsilon;
    a2 = 1. / (2. * epsilon) * a1;
    b2 = a1;
    c2 = 0.5 * a1 * epsilon;
    d2 = 1. / 6. * a1 * epsilon * epsilon;
}

double SoftUnilateralConstraint::evaluate(double x) const {
    x = x - limit;
    if (x > 0)
        return 0.5 * a1 * x * x + b1 * x + c1;
    if (x > -epsilon)
        return 1.0 / 3 * a2 * x * x * x + 0.5 * b2 * x * x + c2 * x + d2;
    return 0;
}

double SoftUnilateralConstraint::computeDerivative(double x) const {
    x = x - limit;
    if (x > 0)
        return a1 * x + b1;
    if (x > -epsilon)
        return a2 * x * x + b2 * x + c2;
    return 0;
}

double SoftUnilateralConstraint::computeSecondDerivative(double x) const {
    x = x - limit;
    if (x > 0)
        return a1;
    if (x > -epsilon)
        return 2 * a2 * x + b2;
    return 0;
}

}  // namespace Opt
}  // namespace DCA
/**
 * This file holds a soft upper limit constraint.
 * 
 * @author Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_SOFTUPPERLIMITCONSTRAINT_H__
#define __DCA_SOFTUPPERLIMITCONSTRAINT_H__

class SoftUpperLimitConstraint {
private:
    double a1, b1, c1, a2, b2, c2, d2, epsilon;
    double limit = 0;

public:
    SoftUpperLimitConstraint(double l, double stiffness, double epsilon) {
        this->limit = l;
        this->epsilon = epsilon;
        a1 = stiffness;

        b1 = 0.5 * a1 * epsilon;
        c1 = 1. / 6. * a1 * epsilon * epsilon;
        a2 = 1. / (2. * epsilon) * a1;
        b2 = a1;
        c2 = 0.5 * a1 * epsilon;
        d2 = 1. / 6. * a1 * epsilon * epsilon;
    }

    virtual ~SoftUpperLimitConstraint() {}

    void setLimit(double l) { limit = l; }
    void setEpsilon(double eps) { epsilon = eps; }

    double compute_F(double x) const {
        x = x - limit;
        if (x > 0) return 0.5 * a1 * x * x + b1 * x + c1;
        if (x > -epsilon)
            return 1.0 / 3 * a2 * x * x * x + 0.5 * b2 * x * x + c2 * x + d2;
        return 0;
    }

    double compute_dFdX(double x) const {
        x = x - limit;
        if (x > 0) return a1 * x + b1;
        if (x > -epsilon) return a2 * x * x + b2 * x + c2;
        return 0;
    }

    double compute_d2FdX2(double x) const {
        x = x - limit;
        if (x > 0) return a1;
        if (x > -epsilon) return 2 * a2 * x + b2;
        return 0;
    }
};

#endif /* __DCA_SOFTUPPERLIMITCONSTRAINT_H__ */
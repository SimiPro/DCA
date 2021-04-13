/**
 * This file holds the implementation for each shape primitive.
 * 
 * @author Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_PRIMITIVES_H__
#define __DCA_PRIMITIVES_H__

#include "utils.h"

namespace DCA {

/**
 * Shape Primitive: Sphere
 */
class Sphere : public CollisionWithAll {
public:
    Sphere(const Vector3d &position, const double &radius)
        : m_position(position), m_radius(radius) {}

    /**
     * Sphere vs. Sphere
     */
    virtual double compute_D(const Sphere &other) const override {
        return (m_position - other.m_position).norm() -
               (m_radius + other.m_radius);
    }

    virtual void compute_dDdP(VectorXd &dDdP,
                              const Sphere &other) const override {
        Vector3d v = m_position - other.m_position;
        double v_norm = v.norm();

        if (v_norm < EPSILON) {
            v_norm = EPSILON;
        }

        Vector3d g = v / v_norm;

        dDdP.resize(6);
        dDdP.head(3) = g;
        dDdP.tail(3) = -g;
    }

    virtual void compute_d2DdP2(MatrixXd &d2DdP2,
                                const Sphere &other) const override {
        Vector3d v = m_position - other.m_position;
        double v_norm = v.norm();

        if (v_norm < EPSILON) {
            v_norm = EPSILON;
        }

        Matrix3d h =
            (Matrix3d::Identity() * v_norm - (v * v.transpose()) / v_norm) /
            (v_norm * v_norm);

        d2DdP2.resize(6, 6);
        d2DdP2.block(0, 0, 3, 3) = h;
        d2DdP2.block(3, 3, 3, 3) = h;
        d2DdP2.block(0, 3, 3, 3) = -h;
        d2DdP2.block(3, 0, 3, 3) = -h;
    }

    // Sphere vs. Capsule
    virtual double compute_D(const Capsule &other) const { return -1.; }
    compute_dDdP_FD(Capsule);
    compute_d2DdP2_FD(Capsule);

    // Helpers
    Vector3d getPosition() const { return m_position; }
    double getRadius() const { return m_radius; }

private:
    Vector3d m_position;
    double m_radius;
};

class Capsule : public CollisionWithAll {
public:
    /**
     * A capsule is defined by two points and a radius.
     */
    Capsule(const Vector3d &startPosition, const Vector3d &endPosition,
            const double &radius)
        : m_startPosition(startPosition),
          m_endPosition(endPosition),
          m_radius(radius) {}

    virtual double compute_D(const Sphere &other) const { return -1.; }
    virtual double compute_D(const Capsule &other) const { return -1.; }

    compute_dDdP_FD(Sphere);
    compute_dDdP_FD(Capsule);

    compute_d2DdP2_FD(Sphere);
    compute_d2DdP2_FD(Capsule);

    // Helpers
    Vector3d getStartPosition() const { return m_startPosition; }
    Vector3d getEndPosition() const { return m_endPosition; }
    double getRadius() const { return m_radius; }
private:
    Vector3d m_startPosition;
    Vector3d m_endPosition;
    double m_radius;
};

// ========== CUSTOM FUNCTIONS ==========
// The following three functions are for easier access from outside
// such that the user does not need to write custom std::visit calls.

/**
 * Compute the distance of a given pair (indices) from the primitives
 */
double compute_D(const pair_t &pair,
                 const std::vector<primitive_t> &primitives) {
    return std::visit(
        [](const auto &cp1, const auto &cp2) { return cp1.compute_D(cp2); },
        primitives.at(pair.first), primitives.at(pair.second));
}

/**
 * Computes the gradient of a given pair and primitives
 */
void compute_dDdP(VectorXd &dDdP, const pair_t &pair,
                  const std::vector<primitive_t> &primitives) {
    std::visit([&dDdP](const auto &cp1,
                       const auto &cp2) { cp1.compute_dDdP(dDdP, cp2); },
               primitives.at(pair.first), primitives.at(pair.second));
}

/**
 * Computes the hessian of a given pair and primitives
 */
void compute_d2DdP2(MatrixXd &d2DdP2, const pair_t &pair,
                    const std::vector<primitive_t> &primitives) {
    std::visit([&d2DdP2](const auto &cp1,
                         const auto &cp2) { cp1.compute_d2DdP2(d2DdP2, cp2); },
               primitives.at(pair.first), primitives.at(pair.second));
}

}  // namespace DCA
#endif /* __DCA_PRIMITIVES_H__ */
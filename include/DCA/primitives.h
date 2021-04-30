/**
 * This file holds the implementation for each shape primitive.
 * 
 * @author Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_PRIMITIVES_H__
#define __DCA_PRIMITIVES_H__

#include <variant>

#include "utils.h"

namespace DCA {

struct Sphere {
    Sphere(const Vector3d& pos, const double& r) : position(pos), radius(r) {}
    Vector3d position;
    double radius;
};

struct Capsule {
    Capsule(const Vector3d& startPos, const Vector3d& endPos, const double& r)
        : startPosition(startPos), endPosition(endPos), radius(r) {}
    Vector3d startPosition;
    Vector3d endPosition;
    double radius;
};

using primitive_t = std::variant<Sphere, Capsule>;

}  // namespace DCA
#endif /* __DCA_PRIMITIVES_H__ */
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
#include <variant>

namespace DCA {

    struct Sphere {
        Vector3d position;
        double radius;
    };

    struct Capsule {
        Vector3d startPosition;
        Vector3d endPosition;
        double radius;
    }


    using primitive_t = std::variant<Sphere, Capsule>;

}  // namespace DCA
#endif /* __DCA_PRIMITIVES_H__ */
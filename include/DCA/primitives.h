#ifndef __DCA_PRIMITIVES_H__
#define __DCA_PRIMITIVES_H__

#include <variant>

#include "utils.h"

namespace DCA {

/**
 * @brief Definition of a %Sphere
 * 
 * A Sphere has a center (position) and a radius.
 */
struct Sphere {
    /**
     * @brief Construct a new %Sphere from a given position and radius.
     * @param[in] pos The position of the sphere.
     * @param[in] r The radius of the sphere.
     */
    Sphere(const Vector3d& pos, const double& r) : position(pos), radius(r) {}
    Vector3d position;  ///< The internal storage for the position.
    double radius;      ///< The internal storage for the radius.
};

/**
 * @brief Definition of a %Capsule
 * 
 * A %Capsule is defined by two points (start- and end-position) and a radius.
 */
struct Capsule {
    /**
     * @brief Construct a new %Capsule from a given positions and radius.
     * @param[in] startPos The start-position of the capsule.
     * @param[in] endPos The end-position of the capsule.
     * @param[in] r The radius of the sphere.
     */
    Capsule(const Vector3d& startPos, const Vector3d& endPos, const double& r)
        : startPosition(startPos), endPosition(endPos), radius(r) {}
    Vector3d startPosition;  ///< The internal storage for the startPosition.
    Vector3d endPosition;    ///< The internal storage for the endPosition.
    double radius;           ///< The internal storage for the radius.
};

struct Plane {
    Plane(const Vector3d& point, const Vector3d& normal)
        : point(point), normal(normal.normalized()) {}
    Plane(const Vector3d& p1, const Vector3d& p2, const Vector3d& p3)
        : point(p1), normal((p2 - p1).cross(p3 - p2).normalized()) {}
    Vector3d point;
    Vector3d normal;
};

/**
 * @brief All possible primitives.
 * This helper holds a variant of all possible primitives supported by this library.
 */
using primitive_t = std::variant<Sphere, Capsule, Plane>;

}  // namespace DCA
#endif /* __DCA_PRIMITIVES_H__ */
#include <DCA/API.h>
#include <DCA/Interactions/CapsuleVsCapsule.h>
#include <DCA/Interactions/PlaneVsCapsule.h>
#include <DCA/Interactions/PlaneVsPlane.h>
#include <DCA/Interactions/PlaneVsSphere.h>
#include <DCA/Interactions/SphereVsCapsule.h>
#include <DCA/Interactions/SphereVsSphere.h>

#define READ_SPHERE_SPHERE_DATA(s_a, s_b) \
    Vector2d props;                       \
    props(0) = s_a.radius;                \
    props(1) = s_b.radius;                \
    Vector6d P;                           \
    P << s_a.position, s_b.position;

#define READ_SPHERE_CAPSULE_DATA(s_a, c_b) \
    Vector2d props;                        \
    props(0) = s_a.radius;                 \
    props(1) = c_b.radius;                 \
    Vector9d P;                            \
    P << s_a.position, c_b.startPosition, c_b.endPosition;

#define READ_CAPSULE_CAPSULE_DATA(c_a, c_b) \
    Vector2d props;                         \
    props(0) = c_a.radius;                  \
    props(1) = c_b.radius;                  \
    Vector12d P;                            \
    P << c_a.startPosition, c_a.endPosition, c_b.startPosition, c_b.endPosition;

#define READ_PLANE_PLANE_DATA(p_a, p_b)                \
    Vector12d P;                                       \
    P << p_a.point, p_a.normal, p_b.point, p_b.normal; \
    Vector0d props;

#define READ_PLANE_CAPSULE_DATA(p_a, c_b) \
    Vector1d props;                       \
    props(0) = c_b.radius;                \
    Vector12d P;                          \
    P << p_a.point, p_a.normal, c_b.startPosition, c_b.endPosition;

#define READ_PLANE_SPHERE_DATA(p_a, s_b) \
    Vector1d props;                      \
    props(0) = s_b.radius;               \
    Vector9d P;                          \
    P << p_a.point, p_a.normal, s_b.position;

namespace DCA {

double API::compute_D(const Sphere& s_a, const Sphere& s_b) {
    READ_SPHERE_SPHERE_DATA(s_a, s_b);
    return SphereVsSphere::compute_D(P, props);
}

void API::compute_dDdP(Vector6d& dDdP, const Sphere& s_a, const Sphere& s_b) {
    READ_SPHERE_SPHERE_DATA(s_a, s_b);
    SphereVsSphere::compute_dDdP(dDdP, P, props);
}

void API::compute_dDdP(Vector3d& dDdP, const Sphere& s_a, const Sphere& s_b) {
    READ_SPHERE_SPHERE_DATA(s_a, s_b);
    SphereVsSphere::compute_dDdP(dDdP, P, props);
}

void API::compute_d2DdP2(Matrix6d& d2DdP2, const Sphere& s_a,
                         const Sphere& s_b) {
    READ_SPHERE_SPHERE_DATA(s_a, s_b);
    SphereVsSphere::compute_d2DdP2(d2DdP2, P, props);
}

void API::compute_d2DdP2(Matrix3d& d2DdP2, const Sphere& s_a,
                         const Sphere& s_b) {
    READ_SPHERE_SPHERE_DATA(s_a, s_b);
    SphereVsSphere::compute_d2DdP2(d2DdP2, P, props);
}

double API::compute_D(const Sphere& s_a, const Capsule& c_b) {
    READ_SPHERE_CAPSULE_DATA(s_a, c_b);
    return SphereVsCapsule::compute_D(P, props);
}

void API::compute_dDdP(Vector9d& dDdP, const Sphere& s_a, const Capsule& c_b) {
    READ_SPHERE_CAPSULE_DATA(s_a, c_b);
    SphereVsCapsule::compute_dDdP(dDdP, P, props);
}

void API::compute_dDdP(Vector3d& dDdP, const Sphere& s_a, const Capsule& c_b) {
    READ_SPHERE_CAPSULE_DATA(s_a, c_b);
    SphereVsCapsule::compute_dDdP(dDdP, P, props);
}

void API::compute_dDdP(Vector6d& dDdP, const Sphere& s_a, const Capsule& c_b) {
    READ_SPHERE_CAPSULE_DATA(s_a, c_b);
    SphereVsCapsule::compute_dDdP(dDdP, P, props);
}

void API::compute_d2DdP2(Matrix9d& d2DdP2, const Sphere& s_a,
                         const Capsule& c_b) {
    READ_SPHERE_CAPSULE_DATA(s_a, c_b);
    SphereVsCapsule::compute_d2DdP2(d2DdP2, P, props);
}

void API::compute_d2DdP2(Matrix3d& d2DdP2, const Sphere& s_a,
                         const Capsule& c_b) {
    READ_SPHERE_CAPSULE_DATA(s_a, c_b);
    SphereVsCapsule::compute_d2DdP2(d2DdP2, P, props);
}

void API::compute_d2DdP2(Matrix6d& d2DdP2, const Sphere& s_a,
                         const Capsule& c_b) {
    READ_SPHERE_CAPSULE_DATA(s_a, c_b);
    SphereVsCapsule::compute_d2DdP2(d2DdP2, P, props);
}

double API::compute_D(const Capsule& c_a, const Sphere& s_b) {
    READ_SPHERE_CAPSULE_DATA(s_b, c_a);
    return SphereVsCapsule::compute_D(P, props);
}

void API::compute_dDdP(Vector9d& dDdP, const Capsule& c_a, const Sphere& s_b) {
    READ_SPHERE_CAPSULE_DATA(s_b, c_a);
    Vector9d dDdP_flip;
    SphereVsCapsule::compute_dDdP(dDdP_flip, P, props);
    dDdP.head(6) = dDdP_flip.tail(6);
    dDdP.tail(3) = dDdP_flip.head(3);
}

void API::compute_dDdP(Vector6d& dDdP, const Capsule& c_a, const Sphere& s_b) {
    READ_SPHERE_CAPSULE_DATA(s_b, c_a);
    SphereVsCapsule::compute_dDdP(dDdP, P, props);
}

void API::compute_dDdP(Vector3d& dDdP, const Capsule& c_a, const Sphere& s_b) {
    READ_SPHERE_CAPSULE_DATA(s_b, c_a);
    SphereVsCapsule::compute_dDdP(dDdP, P, props);
}

void API::compute_d2DdP2(Matrix9d& d2DdP2, const Capsule& c_a,
                         const Sphere& s_b) {
    READ_SPHERE_CAPSULE_DATA(s_b, c_a);
    Matrix9d d2DdP2_flip;
    SphereVsCapsule::compute_d2DdP2(d2DdP2_flip, P, props);
    d2DdP2.block(0, 0, 6, 6) = d2DdP2_flip.block(3, 3, 6, 6);
    d2DdP2.block(6, 6, 3, 3) = d2DdP2_flip.block(0, 0, 3, 3);

    d2DdP2.block(0, 6, 6, 3) = d2DdP2_flip.block(3, 0, 6, 3);
    d2DdP2.block(6, 0, 3, 6) = d2DdP2_flip.block(0, 3, 3, 6);
}

void API::compute_d2DdP2(Matrix6d& d2DdP2, const Capsule& c_a,
                         const Sphere& s_b) {
    READ_SPHERE_CAPSULE_DATA(s_b, c_a);
    SphereVsCapsule::compute_d2DdP2(d2DdP2, P, props);
}

void API::compute_d2DdP2(Matrix3d& d2DdP2, const Capsule& c_a,
                         const Sphere& s_b) {
    READ_SPHERE_CAPSULE_DATA(s_b, c_a);
    SphereVsCapsule::compute_d2DdP2(d2DdP2, P, props);
}

double API::compute_D(const Capsule& c_a, const Capsule& c_b) {
    READ_CAPSULE_CAPSULE_DATA(c_a, c_b);
    return CapsuleVsCapsule::compute_D(P, props);
}

void API::compute_dDdP(Vector12d& dDdP, const Capsule& c_a,
                       const Capsule& c_b) {
    READ_CAPSULE_CAPSULE_DATA(c_a, c_b);
    CapsuleVsCapsule::compute_dDdP(dDdP, P, props);
}

void API::compute_dDdP(Vector6d& dDdP, const Capsule& c_a, const Capsule& c_b) {
    READ_CAPSULE_CAPSULE_DATA(c_a, c_b);
    CapsuleVsCapsule::compute_dDdP(dDdP, P, props);
}

void API::compute_d2DdP2(Matrix12d& d2DdP2, const Capsule& c_a,
                         const Capsule& c_b) {
    READ_CAPSULE_CAPSULE_DATA(c_a, c_b);
    CapsuleVsCapsule::compute_d2DdP2(d2DdP2, P, props);
}

void API::compute_d2DdP2(Matrix6d& d2DdP2, const Capsule& c_a,
                         const Capsule& c_b) {
    READ_CAPSULE_CAPSULE_DATA(c_a, c_b);
    CapsuleVsCapsule::compute_d2DdP2(d2DdP2, P, props);
}

double API::compute_D(const Plane& p_a, const Plane& p_b) {
    READ_PLANE_PLANE_DATA(p_a, p_b);
    return PlaneVsPlane::compute_D(P, props);
}

void API::compute_dDdP(Vector12d& dDdP, const Plane& p_a, const Plane& p_b) {
    READ_PLANE_PLANE_DATA(p_a, p_b);
    PlaneVsPlane::compute_dDdP(dDdP, P, props);
}

void API::compute_d2DdP2(Matrix12d& d2DdP2, const Plane& p_a,
                         const Plane& p_b) {
    READ_PLANE_PLANE_DATA(p_a, p_b);
    PlaneVsPlane::compute_d2DdP2(d2DdP2, P, props);
}

void API::compute_dDdP(Vector6d& dDdP, const Plane& p_a, const Plane& p_b) {
    READ_PLANE_PLANE_DATA(p_a, p_b);
    PlaneVsPlane::compute_dDdP(dDdP, P, props);
}

void API::compute_d2DdP2(Matrix6d& d2DdP2, const Plane& p_a, const Plane& p_b) {
    READ_PLANE_PLANE_DATA(p_a, p_b);
    PlaneVsPlane::compute_d2DdP2(d2DdP2, P, props);
}

double API::compute_D(const Plane& p_a, const Capsule& c_b) {
    READ_PLANE_CAPSULE_DATA(p_a, c_b);
    return PlaneVsCapsule::compute_D(P, props);
}

void API::compute_dDdP(Vector12d& dDdP, const Plane& p_a, const Capsule& c_b) {
    READ_PLANE_CAPSULE_DATA(p_a, c_b);
    PlaneVsCapsule::compute_dDdP(dDdP, P, props);
}

void API::compute_d2DdP2(Matrix12d& d2DdP2, const Plane& p_a,
                         const Capsule& c_b) {
    READ_PLANE_CAPSULE_DATA(p_a, c_b);
    PlaneVsCapsule::compute_d2DdP2(d2DdP2, P, props);
}

void API::compute_dDdP(Vector6d& dDdP, const Plane& p_a, const Capsule& c_b) {
    READ_PLANE_CAPSULE_DATA(p_a, c_b);
    PlaneVsCapsule::compute_dDdP(dDdP, P, props);
}

void API::compute_d2DdP2(Matrix6d& d2DdP2, const Plane& p_a,
                         const Capsule& c_b) {
    READ_PLANE_CAPSULE_DATA(p_a, c_b);
    PlaneVsCapsule::compute_d2DdP2(d2DdP2, P, props);
}

double API::compute_D(const Capsule& c_a, const Plane& p_b) {
    READ_PLANE_CAPSULE_DATA(p_b, c_a);
    return PlaneVsCapsule::compute_D(P, props);
}

void API::compute_dDdP(Vector12d& dDdP, const Capsule& c_a, const Plane& p_b) {
    READ_PLANE_CAPSULE_DATA(p_b, c_a);
    Vector12d dDdP_flipped;
    PlaneVsCapsule::compute_dDdP(dDdP_flipped, P, props);
    dDdP.head(6) = dDdP_flipped.tail(6);
    dDdP.tail(6) = dDdP_flipped.head(6);
}

void API::compute_d2DdP2(Matrix12d& d2DdP2, const Capsule& c_a,
                         const Plane& p_b) {
    READ_PLANE_CAPSULE_DATA(p_b, c_a);
    Matrix12d d2DdP2_flipped;
    PlaneVsCapsule::compute_d2DdP2(d2DdP2_flipped, P, props);
    d2DdP2.block(0, 0, 6, 6) = d2DdP2_flipped.block(6, 6, 6, 6);
    d2DdP2.block(6, 6, 6, 6) = d2DdP2_flipped.block(0, 0, 6, 6);

    d2DdP2.block(0, 6, 6, 6) = d2DdP2_flipped.block(6, 0, 6, 6);
    d2DdP2.block(6, 0, 6, 6) = d2DdP2_flipped.block(0, 6, 6, 6);
}

void API::compute_dDdP(Vector6d& dDdP, const Capsule& c_a, const Plane& p_b) {
    READ_PLANE_CAPSULE_DATA(p_b, c_a);
    Vector12d dDdP_flipped;
    PlaneVsCapsule::compute_dDdP(dDdP_flipped, P, props);
    dDdP = dDdP_flipped.tail(6);
}

void API::compute_d2DdP2(Matrix6d& d2DdP2, const Capsule& c_a,
                         const Plane& p_b) {
    READ_PLANE_CAPSULE_DATA(p_b, c_a);
    Matrix12d d2DdP2_flipped;
    PlaneVsCapsule::compute_d2DdP2(d2DdP2_flipped, P, props);
    d2DdP2 = d2DdP2_flipped.block(6, 6, 6, 6);
}

double API::compute_D(const Plane& p_a, const Sphere& s_b) {
    READ_PLANE_SPHERE_DATA(p_a, s_b);
    return PlaneVsSphere::compute_D(P, props);
}

void API::compute_dDdP(Vector9d& dDdP, const Plane& p_a, const Sphere& s_b) {
    READ_PLANE_SPHERE_DATA(p_a, s_b);
    PlaneVsSphere::compute_dDdP(dDdP, P, props);
}

void API::compute_d2DdP2(Matrix9d& d2DdP2, const Plane& p_a,
                         const Sphere& s_b) {
    READ_PLANE_SPHERE_DATA(p_a, s_b);
    PlaneVsSphere::compute_d2DdP2(d2DdP2, P, props);
}

void API::compute_dDdP(Vector6d& dDdP, const Plane& p_a, const Sphere& s_b) {
    READ_PLANE_SPHERE_DATA(p_a, s_b);
    PlaneVsSphere::compute_dDdP(dDdP, P, props);
}

void API::compute_d2DdP2(Matrix6d& d2DdP2, const Plane& p_a,
                         const Sphere& s_b) {
    READ_PLANE_SPHERE_DATA(p_a, s_b);
    PlaneVsSphere::compute_d2DdP2(d2DdP2, P, props);
}

void API::compute_dDdP(Vector3d& dDdP, const Plane& p_a, const Sphere& s_b) {
    READ_PLANE_SPHERE_DATA(p_a, s_b);
    PlaneVsSphere::compute_dDdP(dDdP, P, props);
}

void API::compute_d2DdP2(Matrix3d& d2DdP2, const Plane& p_a,
                         const Sphere& s_b) {
    READ_PLANE_SPHERE_DATA(p_a, s_b);
    PlaneVsSphere::compute_d2DdP2(d2DdP2, P, props);
}

double API::compute_D(const Sphere& s_a, const Plane& p_b) {
    READ_PLANE_SPHERE_DATA(p_b, s_a);
    return PlaneVsSphere::compute_D(P, props);
}

void API::compute_dDdP(Vector9d& dDdP, const Sphere& s_a, const Plane& p_b) {
    READ_PLANE_SPHERE_DATA(p_b, s_a);
    Vector9d dDdP_flipped;
    PlaneVsSphere::compute_dDdP(dDdP_flipped, P, props);
    dDdP.head(3) = dDdP_flipped.tail(3);
    dDdP.tail(6) = dDdP_flipped.head(6);
}

void API::compute_d2DdP2(Matrix9d& d2DdP2, const Sphere& s_a,
                         const Plane& p_b) {
    READ_PLANE_SPHERE_DATA(p_b, s_a);
    PlaneVsSphere::compute_d2DdP2(d2DdP2, P, props);
    Matrix9d d2DdP2_flipped;
    PlaneVsSphere::compute_d2DdP2(d2DdP2_flipped, P, props);
    d2DdP2.block(0, 0, 3, 3) = d2DdP2_flipped.block(6, 6, 3, 3);
    d2DdP2.block(3, 3, 6, 6) = d2DdP2_flipped.block(0, 0, 6, 6);

    d2DdP2.block(0, 3, 3, 6) = d2DdP2_flipped.block(6, 0, 3, 6);
    d2DdP2.block(3, 0, 6, 3) = d2DdP2_flipped.block(0, 6, 6, 3);
}

void API::compute_dDdP(Vector6d& dDdP, const Sphere& s_a, const Plane& p_b) {
    compute_dDdP(dDdP, p_b, s_a);  // can call this one
}

void API::compute_d2DdP2(Matrix6d& d2DdP2, const Sphere& s_a,
                         const Plane& p_b) {
    compute_d2DdP2(d2DdP2, p_b, s_a);  // can call this one
}

void API::compute_dDdP(Vector3d& dDdP, const Sphere& s_a, const Plane& p_b) {
    compute_dDdP(dDdP, p_b, s_a);  // can call this one
}

void API::compute_d2DdP2(Matrix3d& d2DdP2, const Sphere& s_a,
                         const Plane& p_b) {
    compute_d2DdP2(d2DdP2, p_b, s_a);  // can call this one
}

double API::compute_D(const primitive_t& p_a, const primitive_t& p_b) {
    return std::visit(overloaded{[&](const Sphere& s_a, const Sphere& s_b) {
                                     return compute_D(s_a, s_b);
                                 },
                                 [&](const Sphere& s_a, const Capsule& c_b) {
                                     return compute_D(s_a, c_b);
                                 },
                                 [&](const Capsule& c_a, const Sphere& s_b) {
                                     return compute_D(c_a, s_b);
                                 },
                                 [&](const Capsule& c_a, const Capsule& c_b) {
                                     return compute_D(c_a, c_b);
                                 },
                                 [&](const Plane& p_a, const Plane& p_b) {
                                     return compute_D(p_a, p_b);
                                 },
                                 [&](const Plane& p_a, const Capsule& c_b) {
                                     return compute_D(p_a, c_b);
                                 },
                                 [&](const Capsule& c_a, const Plane& p_b) {
                                     return compute_D(c_a, p_b);
                                 },
                                 [&](const Plane& p_a, const Sphere& s_b) {
                                     return compute_D(p_a, s_b);
                                 },
                                 [&](const Sphere& s_a, const Plane& p_b) {
                                     return compute_D(s_a, p_b);
                                 }},
                      p_a, p_b);
}

void API::compute_dDdP(VectorXd& dDdP, const primitive_t& p_a,
                       const primitive_t& p_b) {
    std::visit(overloaded{[&](const Sphere& s_a, const Sphere& s_b) {
                              Vector6d dDdP_full;
                              compute_dDdP(dDdP_full, s_a, s_b);
                              dDdP = dDdP_full;
                          },
                          [&](const Sphere& s_a, const Capsule& c_b) {
                              Vector9d dDdP_full;
                              compute_dDdP(dDdP_full, s_a, c_b);
                              dDdP = dDdP_full;
                          },
                          [&](const Capsule& c_a, const Sphere& s_b) {
                              Vector9d dDdP_full;
                              compute_dDdP(dDdP_full, c_a, s_b);
                              dDdP = dDdP_full;
                          },
                          [&](const Capsule& c_a, const Capsule& c_b) {
                              Vector12d dDdP_full;
                              compute_dDdP(dDdP_full, c_a, c_b);
                              dDdP = dDdP_full;
                          },
                          [&](const Plane& p_a, const Plane& p_b) {
                              Vector12d dDdP_full;
                              compute_dDdP(dDdP_full, p_a, p_b);
                              dDdP = dDdP_full;
                          },
                          [&](const Plane& p_a, const Capsule& c_b) {
                              Vector12d dDdP_full;
                              compute_dDdP(dDdP_full, p_a, c_b);
                              dDdP = dDdP_full;
                          },
                          [&](const Capsule& c_a, const Plane& p_b) {
                              Vector12d dDdP_full;
                              compute_dDdP(dDdP_full, c_a, p_b);
                              dDdP = dDdP_full;
                          },
                          [&](const Plane& p_a, const Sphere& s_b) {
                              Vector9d dDdP_full;
                              compute_dDdP(dDdP_full, p_a, s_b);
                              dDdP = dDdP_full;
                          },
                          [&](const Sphere& s_a, const Plane& p_b) {
                              Vector9d dDdP_full;
                              compute_dDdP(dDdP_full, s_a, p_b);
                              dDdP = dDdP_full;
                          }},
               p_a, p_b);
}

void API::compute_d2DdP2(MatrixXd& d2DdP2, const primitive_t& p_a,
                         const primitive_t& p_b) {
    std::visit(overloaded{[&](const Sphere& s_a, const Sphere& s_b) {
                              Matrix6d d2DdP2_full;
                              compute_d2DdP2(d2DdP2_full, s_a, s_b);
                              d2DdP2 = d2DdP2_full;
                          },
                          [&](const Sphere& s_a, const Capsule& c_b) {
                              Matrix9d d2DdP2_full;
                              compute_d2DdP2(d2DdP2_full, s_a, c_b);
                              d2DdP2 = d2DdP2_full;
                          },
                          [&](const Capsule& c_a, const Sphere& s_b) {
                              Matrix9d d2DdP2_full;
                              compute_d2DdP2(d2DdP2_full, c_a, s_b);
                              d2DdP2 = d2DdP2_full;
                          },
                          [&](const Capsule& c_a, const Capsule& c_b) {
                              Matrix12d d2DdP2_full;
                              compute_d2DdP2(d2DdP2_full, c_a, c_b);
                              d2DdP2 = d2DdP2_full;
                          },
                          [&](const Plane& p_a, const Plane& p_b) {
                              Matrix12d d2DdP2_full;
                              compute_d2DdP2(d2DdP2_full, p_a, p_b);
                              d2DdP2 = d2DdP2_full;
                          },
                          [&](const Plane& p_a, const Capsule& c_b) {
                              Matrix12d d2DdP2_full;
                              compute_d2DdP2(d2DdP2_full, p_a, c_b);
                              d2DdP2 = d2DdP2_full;
                          },
                          [&](const Capsule& c_a, const Plane& p_b) {
                              Matrix12d d2DdP2_full;
                              compute_d2DdP2(d2DdP2_full, c_a, p_b);
                              d2DdP2 = d2DdP2_full;
                          },
                          [&](const Plane& p_a, const Sphere& s_b) {
                              Matrix9d d2DdP2_full;
                              compute_d2DdP2(d2DdP2_full, p_a, s_b);
                              d2DdP2 = d2DdP2_full;
                          },
                          [&](const Sphere& s_a, const Plane& p_b) {
                              Matrix9d d2DdP2_full;
                              compute_d2DdP2(d2DdP2_full, s_a, p_b);
                              d2DdP2 = d2DdP2_full;
                          }},
               p_a, p_b);
}

double API::compute_D(const pair_t& pair,
                      const std::vector<primitive_t>& primitives) {
    return compute_D(primitives.at(pair.first), primitives.at(pair.second));
}

void API::compute_dDdP(VectorXd& dDdP, const pair_t& pair,
                       const std::vector<primitive_t>& primitives) {
    compute_dDdP(dDdP, primitives.at(pair.first), primitives.at(pair.second));
}

void API::compute_d2DdP2(MatrixXd& d2DdP2, const pair_t& pair,
                         const std::vector<primitive_t>& primitives) {
    compute_d2DdP2(d2DdP2, primitives.at(pair.first),
                   primitives.at(pair.second));
}

}  // namespace DCA
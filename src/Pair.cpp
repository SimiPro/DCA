#include <DCA/Pair.h>
#include <DCA/Interactions/PlaneVsSphere.h>

namespace DCA {

std::vector<pair_t> PermutationPairGenerator::generate(
    const std::vector<primitive_t> &primitives) const {
    // Create all possible permutations.
    std::vector<pair_t> ret;
    ret.reserve(primitives.size() * primitives.size() - primitives.size());

    for (size_t i = 0; i < primitives.size(); i++) {
        for (size_t j = i + 1; j < primitives.size(); j++) {
            // skip primitve self-collision
            ret.push_back({i, j});
        }
    }

    return ret;
}

NeighborsPairGenerator::NeighborsPairGenerator(const double &radius)
    : m_radius2(radius * radius) {}

std::vector<pair_t> NeighborsPairGenerator::generate(
    const std::vector<primitive_t> &primitives) const {
    std::vector<pair_t> ret;

    std::vector<Vector3d> positions(primitives.size());

    for (size_t i = 0; i < primitives.size(); i++) {
        Vector3d pos = std::visit(
            overloaded{
                [](const Sphere &cp) { return cp.position; },
                [](const Capsule &cp) {
                    return Vector3d(0.5 * (cp.startPosition + cp.endPosition));
                },
                [](const Plane &cp) {
                    return cp.point;  // palceholder. Will be checked later.
                }},
            primitives[i]);
        positions[i] = {pos.x(), pos.y(), pos.z()};
    }

    // now search all pairs where they are in a certain distance
    for (size_t i = 0; i < primitives.size(); i++) {
        // skip primitve self-collision, therefore start at i + 1
        for (size_t j = i + 1; j < primitives.size(); j++) {
            std::visit(
                overloaded{
                    [&](const Plane &p, const Plane &other) {
                        // compare normals, if dot product = 1 or -1 --> dist.
                        // else: 0
                        double normals = p.normal.dot(other.normal);
                        if (fabs(fabs(normals) - 1.) < 1e-8) {
                            // check if distance fullfilled
                            Vector6d P_plane;
                            P_plane << p.point, p.normal;
                            Vector3d pt_on_plane =
                                PlaneVsSphere::getProjectionOfPoint(
                                    P_plane, other.point);
                            if ((pt_on_plane - positions[j]).squaredNorm() <
                                m_radius2) {
                                ret.push_back({i, j});
                            }
                        } else {
                            // intersecting, return pair
                            ret.push_back({i, j});
                        }
                    },
                    [&](const Plane &p, const primitive_t &other) {
                        // map pos. of other to the plane
                        Vector6d P_plane;
                        P_plane << p.point, p.normal;
                        Vector3d pt_on_plane =
                            PlaneVsSphere::getProjectionOfPoint(P_plane,
                                                                positions[j]);
                        // now check distance
                        if ((pt_on_plane - positions[j]).squaredNorm() <
                            m_radius2) {
                            ret.push_back({i, j});
                        }
                    },
                    [&](const primitive_t &a, const primitive_t &b) {
                        // if two "normal" primitives, compare positions
                        if ((positions[i] - positions[j]).squaredNorm() <
                            m_radius2) {
                            ret.push_back({i, j});
                        }
                    }},
                primitives.at(i), primitives.at(j));
        }
    }

    return ret;
}
}  // namespace DCA
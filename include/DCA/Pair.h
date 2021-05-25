#ifndef __DCA_PAIR_H__
#define __DCA_PAIR_H__

#include <DCA/Utils/Primitives.h>
#include <DCA/Utils/Utils.h>

namespace DCA {
/**
 * @brief This is the base class for all generators.
 */
class PairGenerator {
public:
    /**
     * @brief This function generates the pairs given the primitives.
     * 
     * The returned vector consists of pairs, where each pair holds two numbers:
     * The indices of the corresponding primitives which were given.
     * @param[in] primitives All primitives to generate the pairs from.
     * @return A vector of pairs of indices, where each index corresponds to a primitive in the primitives vector.
     */
    virtual std::vector<pair_t> generate(
        const std::vector<primitive_t> &primitives) const = 0;
};

/**
 * @brief This generator creates all possible permutations of pairs.
 * This means, the amount of pairs created is \f$n^2/2\f$, where n is the number of primitives.
 * Pairs are not returned twice (0, 1) and (1, 0).
 */
class PermutationPairGenerator : public PairGenerator {
public:
    /**
     * @copydoc PairGenerator::generate
     */
    virtual std::vector<pair_t> generate(
        const std::vector<primitive_t> &primitives) const override {
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
};

/**
 * This generator computes the pairs which are in a certain threshold from each other.
 * It does so by computing a single position for each primitive and selecting
 * pairs based on the distance.
 */
class NeighborsPairGenerator : public PairGenerator {
public:
    /**
     * @brief Construct this generator with a given radius
     * @param[in] radius The radius to search other primitives in.
     */
    NeighborsPairGenerator(const double &radius) : m_radius2(radius * radius) {}

    /**
     * @copydoc PairGenerator::generate
     */
    virtual std::vector<pair_t> generate(
        const std::vector<primitive_t> &primitives) const override {
        std::vector<pair_t> ret;

        std::vector<Vector3d> positions(primitives.size());

        for (size_t i = 0; i < primitives.size(); i++) {
            Vector3d pos = std::visit(
                overloaded{
                    [](const Sphere &cp) { return cp.position; },
                    [](const Capsule &cp) {
                        return Vector3d(0.5 *
                                        (cp.startPosition + cp.endPosition));
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
                                PlaneVsSphere::getProjectionOfPoint(
                                    P_plane, positions[j]);
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
    }  // namespace DCA

private:
    double m_radius2;  ///< The squared radius
};

}  // namespace DCA
#endif /* __DCA_PAIR_H__ */
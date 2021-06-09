#include <DCA/Interactions/PlaneVsSphere.h>
#include <DCA/Pair.h>
#include <DCA/API.h>

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
    : m_radius(radius) {}

std::vector<pair_t> NeighborsPairGenerator::generate(
    const std::vector<primitive_t> &primitives) const {
    std::vector<pair_t> ret;

    // now search all pairs where they are in a certain distance
    for (size_t i = 0; i < primitives.size(); i++) {
        // skip primitive self-collision, therefore start at i + 1
        for (size_t j = i + 1; j < primitives.size(); j++) {
            double distance = API::compute_D(primitives.at(i), primitives.at(j));
            if (distance < m_radius) {
                ret.push_back({i, j});
            }
        }
    }

    return ret;
}
}  // namespace DCA
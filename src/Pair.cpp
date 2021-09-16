#include <DCA/Pair.h>

namespace DCA {

std::vector<pair_t> PermutationPairGenerator::generate(const std::vector<primitive_t> &primitives) const {
    // Create all possible permutations.
    std::vector<pair_t> ret;
    ret.reserve(primitives.size() * primitives.size() - primitives.size());

    for (size_t i = 0; i < primitives.size(); i++) {
        for (size_t j = i + 1; j < primitives.size(); j++) {
            // skip primitive self-collision
            ret.push_back({i, j});
        }
    }

    return ret;
}

NeighborsPairGenerator::NeighborsPairGenerator(const double &radius) : m_radius(radius) {}

std::vector<pair_t> NeighborsPairGenerator::generate(const std::vector<primitive_t> &primitives) const {
    std::vector<pair_t> ret;

    // now search all pairs where they are in a certain distance
    for (size_t i = 0; i < primitives.size(); i++) {
        // skip primitive self-collision, therefore start at i + 1
        for (size_t j = i + 1; j < primitives.size(); j++) {
            //double distance = API::compute_D(primitives.at(i), primitives.at(j));
            double distance = estimate_distance(primitives.at(i), primitives.at(j));
            if (distance < m_radius) {
                ret.push_back({i, j});
            }
        }
    }

    return ret;
}

std::vector<pair_t> NeighborsPairGenerator::generate(const std::vector<primitive_t> &primitives_a, const std::vector<primitive_t> &primitives_b) const {
    std::vector<pair_t> ret;

    for (size_t i = 0; i < primitives_a.size(); i++) {
        for (size_t j = 0; j < primitives_b.size(); j++) {
            //double distance = API::compute_D(primitives_a.at(i), primitives_b.at(j));
            double distance = estimate_distance(primitives_a.at(i), primitives_b.at(j));
            if (distance < m_radius) {
                ret.push_back({i, j});
            }
        }
    }

    return ret;
}

double NeighborsPairGenerator::estimate_distance(const primitive_t &p_A, const primitive_t &p_B) const {
    Vector3d C_A = std::visit([&](const auto &primitive) { return primitive.get_center_point(); }, p_A);
    Vector3d C_B = std::visit([&](const auto &primitive) { return primitive.get_center_point(); }, p_B);
    double d_A = std::visit([&](const auto &primitive) { return primitive.get_largest_dimension_from_center(); }, p_A);
    double d_B = std::visit([&](const auto &primitive) { return primitive.get_largest_dimension_from_center(); }, p_B);
    return (C_A - C_B).norm() - (d_A + d_B);
};

}  // namespace DCA
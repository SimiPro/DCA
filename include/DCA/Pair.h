#pragma once

#include <DCA/Utils/Primitives.h>

namespace DCA {

/**
 * @brief This generator creates all possible permutations of pairs.
 * This means, the amount of pairs created is \f$n^2/2\f$, where n is the number of primitives.
 * Pairs are not returned twice (0, 1) and (1, 0).
 */
class PermutationPairGenerator {
public:
    /**
     * @brief This function generates all pairs given the primitives.
     * 
     * The returned vector consists of pairs, where each pair holds two numbers:
     * The indices of the corresponding primitives which were given.
     * @param[in] primitives All primitives to generate the pairs from.
     * @return A vector of pairs of indices, where each index corresponds to a primitive in the primitives vector.
     */
    std::vector<pair_t> generate(const std::vector<primitive_t> &primitives) const;
};

/**
 * @brief This generator computes the pairs which are in a certain threshold from each other.
 * It does so by computing a single position for each primitive and selecting
 * pairs based on the distance.
 */
class NeighborsPairGenerator {
public:
    /**
     * @brief Construct this generator with a given radius
     * @param[in] radius The radius to search other primitives in.
     */
    NeighborsPairGenerator(const double &radius);

    /**
     * @brief This function generates the pairs given the primitives, where the pairs are in a certain distance from each other.
     * 
     * The returned vector consists of pairs, where each pair holds two numbers:
     * The indices of the corresponding primitives which were given.
     * @param[in] primitives All primitives to generate the pairs from.
     * @return A vector of pairs of indices, where each index corresponds to a primitive in the primitives vector.
     */
    std::vector<pair_t> generate(const std::vector<primitive_t> &primitives) const;  // namespace DCA

    /**
     * @brief This function generates pairs given two vectors of primitives, where the pairs are in a certain distance from each other and are in separate vectors.
     * 
     * The returned vector consists of pairs, where each pair holds two numbers:
     * The indices of the corresponding primitives which were given. The first index corresponds to the primitives_a vector, the second to the primitives_b vector.
     * @param[in] primitives_a The first set of primitives.
     * @param[in] primitives_b The second set of primitives.
     * @return A vector of pairs of indices, where the first index corresponds to a primitive in the primitives_a vector and the second index to the primitives_b vector.
     * 
     * @attention Not all pairs are reported, only those between primitives_a and primitives_b!
     */
    std::vector<pair_t> generate(const std::vector<primitive_t> &primitives_a, const std::vector<primitive_t> &primitives_b) const;

private:
    /**
    *   @brief Helper function for generate
    */
    double estimate_distance(const primitive_t &p_A, const primitive_t &p_B) const;

private:
    double m_radius;  ///< The radius
};

}  // namespace DCA
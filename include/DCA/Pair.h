#pragma once

#include <DCA/Utils/Primitives.h>

namespace DCA {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
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
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * @brief This generator creates all possible permutations of pairs.
 * This means, the amount of pairs created is \f$n^2/2\f$, where n is the number of primitives.
 * Pairs are not returned twice (0, 1) and (1, 0).
 */
class PermutationPairGenerator : public PairGenerator {
public:
    /**
     * @brief This function generates all pairs given the primitives.
     * 
     * The returned vector consists of pairs, where each pair holds two numbers:
     * The indices of the corresponding primitives which were given.
     * @param[in] primitives All primitives to generate the pairs from.
     * @return A vector of pairs of indices, where each index corresponds to a primitive in the primitives vector.
     */
    virtual std::vector<pair_t> generate(
        const std::vector<primitive_t> &primitives) const override;
};

/**
 * @brief This generator computes the pairs which are in a certain threshold from each other.
 * It does so by computing a single position for each primitive and selecting
 * pairs based on the distance.
 */
class NeighborsPairGenerator : public PairGenerator {
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
    virtual std::vector<pair_t> generate(
        const std::vector<primitive_t> &primitives)
        const override;  // namespace DCA

private:
    double m_radius;  ///< The radius
};

}  // namespace DCA
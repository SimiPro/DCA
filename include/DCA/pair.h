/**
 * This file manages all pairs of collisions.
 * This means, it performs various approaches
 * for a broad collision detection.
 * 
 * @author: Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_PAIR_H__
#define __DCA_PAIR_H__

#include "utils.h"

namespace DCA {
/**
 * This is the base class for all generators.
 */
class PairGenerator {
public:
    /**
     * This function generates the pairs given the primitives.
     * The returned vector consists of pairs, where each pair holds two numbers:
     * The indices of the corresponding primitives which were given.
     */
    virtual std::vector<pair_t> generate(
        const std::vector<primitive_t> &primitives) const = 0;
};

/**
 * This generator creates all possible permutations of pairs.
 * This means, the amount of pairs created is n^2, where n = #primitives.
 */
class PermutationPairGenerator : public PairGenerator {
public:
    virtual std::vector<pair_t> generate(
        const std::vector<primitive_t> &primitives) const override {
        // Create all possible permutations.
        std::vector<pair_t> ret;
        ret.reserve(primitives.size() * primitives.size() - primitives.size());

        for (size_t i = 0; i < primitives.size(); i++) {
            for (size_t j = 0; j < primitives.size(); j++) {
                // skip primitve self-collision
                if (i == j) continue;
                ret.push_back({i, j});
            }
        }

        return ret;
    }
};
}  // namespace DCA
#endif /* __DCA_PAIR_H__ */
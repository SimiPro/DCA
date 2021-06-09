/**
 * This file should serve as an example on how to use this library.
 * 
 * @author: Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#include <iostream>  // for std::cout and std::endl

// After linking the library, it can be included with the following lines:
#include <DCA/API.h>
#include <DCA/Pair.h>

int main(int argc, char const* argv[]) {
    using namespace DCA;

    // A vector of some primitives (using primitive_t)
    std::vector<primitive_t> primitives;

    // Push back some primitives
    primitives.push_back(Sphere{Vector3d(0, 0, 0), 0.5});
    primitives.push_back(Sphere{Vector3d(0, 2, 0), 0.5});
    primitives.push_back(Capsule{Vector3d(1, 0, 0), Vector3d(2, 0, 0), 0.5});

    // Create all possible permutations for the given primitives
    PermutationPairGenerator ppg;
    auto pairs_permutatation = ppg.generate(primitives);
    std::cout << "Permutation reported " << pairs_permutatation.size()
              << " pairs." << std::endl;
    for(auto& p : pairs_permutatation) {
            std::cout << p.first << "/" << p.second << std::endl;
    }

    NeighborsPairGenerator npg(2.);
    auto pairs_neighbors = npg.generate(primitives);

    std::cout << "Neighbors reported " << pairs_neighbors.size() << " pairs."
              << std::endl;
    for(auto& p : pairs_neighbors) {
            std::cout << p.first << "/" << p.second << std::endl;
    }

    // Temporary storage for the gradient and hessian
    VectorXd grad;
    MatrixXd hess;

    // Iterate over each pair and print out the distance and derivatives.
    for (auto& pair : pairs_permutatation) {
        std::cout << "Distance: " << API::compute_D(pair, primitives)
                  << std::endl;

        API::compute_dDdP(grad, pair, primitives);
        std::cout << "Gradient^T:" << std::endl
                  << grad.transpose() << std::endl;

        API::compute_d2DdP2(hess, primitives.at(pair.first),
                            primitives.at(pair.second));
        std::cout << "Hessian:" << std::endl << hess << std::endl;
    }

    return 0;
}

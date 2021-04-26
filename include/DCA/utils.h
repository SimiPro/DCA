/**
 * This file holds utils used all over this library.
 * 
 * @author Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_UTILS_H__
#define __DCA_UTILS_H__

#include <Eigen/Core>
#include <vector>

namespace DCA {

using Eigen::Vector2d;

using Eigen::Matrix3d;
using Eigen::Vector3d;

using Vector6d = Eigen::Matrix<double, 6, 1>;
using Matrix6d = Eigen::Matrix<double, 6, 6>;

using Vector9d = Eigen::Matrix<double, 9, 1>;
using Matrix9d = Eigen::Matrix<double, 9, 9>;

// Easier access to any dynamic matrix
using Eigen::MatrixXd;

// Easier access to any dynamic vector
using Eigen::VectorXd;

// Easier access to a pair (corresponding of two indices)
using pair_t = std::pair<size_t, size_t>;

#ifndef DOXYGEN_SHOULD_SKIP_THIS

// Helper for std::visit
template <class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;
#endif

#define EPSILON 1e-8
#define deltaFD 1e-6

inline double sigmoid(double x, double scale, double shift = 0.5) {
    return 1.0 / (1.0 + std::exp(-1.0 * scale * (x - shift)));
}

}  // namespace DCA

#endif /* __DCA_UTILS_H__ */
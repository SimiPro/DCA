#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

namespace DCA {

typedef unsigned int uint;  ///< Easier access to a unsigned int.

using Vector0d = Eigen::Matrix<double, 0, 1>;    ///< Easier access to a 0 vector of type double.
using Vector1d = Eigen::Matrix<double, 1, 1>;    ///< Easier access to a 1 vector of type double.
using Eigen::Vector2d;                           ///< Easier access to a 2 vector of type double.
using Eigen::Vector3d;                           ///< Easier access to a 3 vector of type double.
using Vector6d = Eigen::Matrix<double, 6, 1>;    ///< Easier access to a 6 vector of type double.
using Vector12d = Eigen::Matrix<double, 12, 1>;  ///< Easier access to a 12 vector.
using Eigen::VectorXd;                           ///< Easier access to any dynamic vector

using Eigen::Matrix3d;                            ///< Easier access to a 3 by 3 matrix.
using Matrix12d = Eigen::Matrix<double, 12, 12>;  ///< Easier access to a 12 by 12 matrix.
using Eigen::MatrixXd;                            ///< Easier access to any dynamic matrix

using pair_t = std::pair<size_t, size_t>;  ///< Easier access to a pair (corresponding of two indices)

// Helper for std::visit, see https://en.cppreference.com/w/cpp/utility/variant/visit
template <class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

}  // namespace DCA

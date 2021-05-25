#ifndef __DCA_UTILS_H__
#define __DCA_UTILS_H__

#include <Eigen/Core>
#include <vector>

namespace DCA {

using Vector0d = Eigen::Matrix<double, 0, 1>; ///< Easier access to a 0 vector.
using Vector1d = Eigen::Matrix<double, 1, 1>; ///< Easier access to a 1 vector.

using Eigen::Vector2d; ///< Easier access to a 2 vector.
using Eigen::Matrix2d; ///< Easier access to a 2 by 2 matrix.

using Eigen::Vector3d; ///< Easier access to a 3 vector.
using Eigen::Matrix3d; ///< Easier access to a 3 by 3 matrix.

using Vector6d = Eigen::Matrix<double, 6, 1>; ///< Easier access to a 6 vector.
using Matrix6d = Eigen::Matrix<double, 6, 6>; ///< Easier access to a 6 by 6 matrix.

using Vector9d = Eigen::Matrix<double, 9, 1>; ///< Easier access to a 9 vector.
using Matrix9d = Eigen::Matrix<double, 9, 9>; ///< Easier access to a 9 by 9 matrix.

using Vector12d = Eigen::Matrix<double, 12, 1>; ///< Easier access to a 12 vector.
using Matrix12d = Eigen::Matrix<double, 12, 12>; ///< Easier access to a 12 by 12 matrix.

using Eigen::VectorXd; ///< Easier access to any dynamic vector
using Eigen::MatrixXd; ///< Easier access to any dynamic matrix

using pair_t = std::pair<size_t, size_t>; ///< Easier access to a pair (corresponding of two indices)

#ifndef DOXYGEN_SHOULD_SKIP_THIS

// Helper for std::visit
template <class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

inline double sigmoid(double x, double scale, double shift = 0.5) {
    return 1.0 / (1.0 + std::exp(-1.0 * scale * (x - shift)));
}

}  // namespace DCA

#endif /* __DCA_UTILS_H__ */
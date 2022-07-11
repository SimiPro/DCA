#pragma once

#include <DCA/Logger.h>
#include <DCA/Utils.h>
#include <DCA/ddR.h>

namespace DCA {

/**
 * @brief Represents Rotations using exponential coordinates.
 * 
 * See e.g. https://arxiv.org/pdf/1312.0788.pdf
 */
class ExpCoords {
public:
    /**
     * @brief Get \f$ \theta \f$ from a rotation matrix.
     * 
     * Returns the exponential map parameterization
     * @param[in] R The rotation matrix
     * @return The exponential map parameterization.
     */
    static Vector3d get_theta(const Matrix3d &R) {
        Vector3d theta;
        theta[0] = R(2, 1) - R(1, 2);
        theta[1] = R(0, 2) - R(2, 0);
        theta[2] = R(1, 0) - R(0, 1);

        if (theta.norm() > 1e-6)
            theta /= theta.norm();

        theta *= safeACOS((R.trace() - 1.0) / 2.0);
        return theta;
    }

    /**
     * @brief Get \f$ \frac{d\theta}{dR_i} \f$.
     * 
     * Returns the derivative of exponential coorinates with respect to a column of R
     * @param[in] R The rotation matrix
     * @param[in] index The column to use. Must be either 0, 1 or 2.
     * @return The exponential map parameterization.
     * 
     * @throws std::runtime_error If the index is out of bounds.
     */
    static Matrix3d get_dThetadRi(const Matrix3d &R, const int &index) {
        if (index < 0 || index > 2)
            throw std::runtime_error("ExpCoords::get_dThetadRi -> invalid index");

        Vector3d theta;
        theta[0] = R(2, 1) - R(1, 2);
        theta[1] = R(0, 2) - R(2, 0);
        theta[2] = R(1, 0) - R(0, 1);

        if (theta.norm() < 1e-6)
            return Matrix3d::Identity();  //This seems to be a good idea...

        const double theta_norm = theta.norm();
        const double theta_norm_3 = theta_norm * theta_norm * theta_norm;
        const double x = (R.trace() - 1.0) / 2.0;
        const double acos = safeACOS(x);
        const double dacos = safeDACOS(x);

        Matrix3d dTdRi;
        if (index == 0) {
            dTdRi.col(0) = 0.5 * dacos * theta / theta_norm;
            dTdRi.col(1) = -acos * (theta * theta[2] / theta_norm_3 - Vector3d::UnitZ() / theta_norm);
            dTdRi.col(2) = acos * (theta * theta[1] / theta_norm_3 - Vector3d::UnitY() / theta_norm);
        } else if (index == 1) {
            dTdRi.col(0) = acos * (theta * theta[2] / theta_norm_3 - Vector3d::UnitZ() / theta_norm);
            dTdRi.col(1) = 0.5 * dacos * theta / theta_norm;
            dTdRi.col(2) = -acos * (theta * theta[0] / theta_norm_3 - Vector3d::UnitX() / theta_norm);
        } else if (index == 2) {
            dTdRi.col(0) = -acos * (theta * theta[1] / theta_norm_3 - Vector3d::UnitY() / theta_norm);
            dTdRi.col(1) = acos * (theta * theta[0] / theta_norm_3 - Vector3d::UnitX() / theta_norm);
            dTdRi.col(2) = 0.5 * dacos * theta / theta_norm;
        }
        return dTdRi;
    }

    /**
     * @brief Get a rotation matrix from given exponential map.
     * @param[in] theta \f$ \theta \f$
     * @return The rotation matrix which is represented by \f$ \theta \f$.
     */
    static Matrix3d get_R(const Vector3d &theta) {
        double t = theta.norm();
        Matrix3d R;
        R(0, 0) = 1;
        R(1, 1) = 1;
        R(2, 2) = 1;
        if (t > tol) {
            double c = (1. - std::cos(t)) / t / t;
            double s = std::sin(t) / t;

            R(0, 0) += c * (theta(0) * theta(0) - t * t);
            R(1, 1) += c * (theta(1) * theta(1) - t * t);
            R(2, 2) += c * (theta(2) * theta(2) - t * t);

            R(0, 1) = -s * theta(2) + c * theta(0) * theta(1);
            R(0, 2) = s * theta(1) + c * theta(0) * theta(2);

            R(1, 0) = s * theta(2) + c * theta(1) * theta(0);
            R(1, 2) = -s * theta(0) + c * theta(1) * theta(2);

            R(2, 0) = -s * theta(1) + c * theta(2) * theta(0);
            R(2, 1) = s * theta(0) + c * theta(2) * theta(1);
            return R;

        } else {
            auto V = getSkewSymmetricMatrix(theta);
            return Matrix3d::Identity() + V + 0.5 * V * V;
        }
    }

    /**
     * @brief Returns the Jacobian of R.
     * 
     * Returns a Jacobian that tells us how the rotation matrix changes with respect to the exponential map parameters that define it.
     * @param[in] theta \f$ \theta \f$
     * @param[in] i The index (column).
     * @return \f$ \frac{dR}{d\theta_i} \f$
     */
    static Matrix3d get_dRi(const Vector3d &theta, int i) {
        Vector3d v = theta;
        Matrix3d V = getSkewSymmetricMatrix(theta);
        double t = theta.norm();
        Matrix3d dR_i;
        Vector3d e_i = Vector3d::Unit(i);
        Matrix3d R = get_R(theta);
        if (t > tol) {
            //derivation here: https://arxiv.org/pdf/1312.0788.pdf - A compact formula for the derivative of the rotation wrt the exponential map parameters that define it
            dR_i = getSkewSymmetricMatrix(Vector3d(v[i] * v + V * (Matrix3d::Identity() - R) * e_i)) / (t * t) * R;
        } else {
            dR_i = getSkewSymmetricMatrix(e_i) + 0.5 * getSkewSymmetricMatrix(e_i) * V + 0.5 * (getSkewSymmetricMatrix(e_i) * V).transpose();
        }
        return dR_i;
    }

    /**
     * @brief Returns the derivative of the Jacobian of R at an index.
     * 
     * Returns the derivative of the Jacobian \f$ \frac{dR}{d\theta_i} \f$ with respect to the exponential map parameters at index j.
     * @param[in] theta \f$ \theta \f$
     * @param[in] i The first index.
     * @param[in] j The second index.
     * @return \f$ \frac{d^2R}{d\theta_i d\theta_j} \f$
     * 
     * @throws std::runtime_error (only release) if either i or j is not in [0,2].
     */
    static Matrix3d get_ddR_i_j(const Vector3d &theta, int i, int j) {
        assert(i >= 0 && i < 3);
        assert(j >= 0 && j < 3);
        Vector3d e_i = Vector3d::Unit(i);
        Vector3d e_j = Vector3d::Unit(j);

        double t = theta.norm();
        if (t > tol) {
            //auto diff'ed code
            if (i == 0 && j == 0)
                return ddR_0_0(theta);
            if (i == 0 && j == 1)
                return ddR_0_1(theta);
            if (i == 0 && j == 2)
                return ddR_0_2(theta);
            if (i == 1 && j == 0)
                return ddR_1_0(theta);
            if (i == 1 && j == 1)
                return ddR_1_1(theta);
            if (i == 1 && j == 2)
                return ddR_1_2(theta);
            if (i == 2 && j == 0)
                return ddR_2_0(theta);
            if (i == 2 && j == 1)
                return ddR_2_1(theta);
            if (i == 2 && j == 2)
                return ddR_2_2(theta);

            throw std::runtime_error("get_ddR_i_j: i and j must be in [0,2]");
        } else {
            return 0.5 * getSkewSymmetricMatrix(e_i) * getSkewSymmetricMatrix(e_j) +
                   0.5 * (getSkewSymmetricMatrix(e_i) * getSkewSymmetricMatrix(e_j)).transpose();
        }

        return Matrix3d::Zero();
    }

    /**
     * @brief Returns world coordinates for a vector.
     * 
     * Returns the world coordinates for the vector x that is expressed in local coordinates.
     * @param[in] theta The exponential map representation.
     * @param[in] x The vector in local coordinates.
     * @return The world coordinates for x.
     */
    static Vector3d get_w(const Vector3d &theta, const Vector3d &x) {
        Vector3d res = get_R(theta) * x;
        return res;
    }

    /**
     * @brief Returns the first derivative of ExpCoords::get_w.
     * 
     * Returns a matrix that tells us how w changes with respect to the orientation.
     * @param[in] theta The exponential map representation.
     * @param[in] x The vector in local coordinates.
     * @return The derivative of w with respect to the orientation.
     */
    static Matrix3d get_dwdr(const Vector3d &theta, const Vector3d &x) {
        Matrix3d dw_dr;
        for (int i = 0; i < 3; ++i)
            dw_dr.col(i) = get_dRi(theta, i) * x;
        return dw_dr;
    }

    /**
     * @brief Returns \f$ \frac{d^2 w}{dR dR_1} \f$
     * 
     * Returns the second derivative of ExpCoords::get_w with respect to the first component of the orientation.
     * @param[in] theta The exponential map representation.
     * @param[in] x The local coordinates of the point.
     * @return The second derivative.
     */
    static Matrix3d get_ddwdr_dr1(const Vector3d &theta, const Vector3d &x) {
        Matrix3d ddwdr_dr1;
        for (int i = 0; i < 3; ++i)
            ddwdr_dr1.col(i) = get_ddR_i_j(theta, i, 0) * x;
        return ddwdr_dr1;
    }

    /**
     * @brief Returns \f$ \frac{d^2 w}{dR dR_2} \f$
     * 
     * Returns the second derivative of ExpCoords::get_w with respect to the second component of the orientation.
     * @param[in] theta The exponential map representation.
     * @param[in] x The local coordinates of the point.
     * @return The second derivative.
     */
    static Matrix3d get_ddwdr_dr2(const Vector3d &theta, const Vector3d &x) {
        Matrix3d ddwdr_dr2;
        for (int i = 0; i < 3; ++i)
            ddwdr_dr2.col(i) = get_ddR_i_j(theta, i, 1) * x;
        return ddwdr_dr2;
    }

    /**
     * @brief Returns \f$ \frac{d^2 w}{dR dR_3} \f$
     * 
     * Returns the second derivative of ExpCoords::get_w with respect to the third component of the orientation.
     * @param[in] theta The exponential map representation.
     * @param[in] x The local coordinates of the point.
     * @return The second derivative.
     */
    static Matrix3d get_ddwdr_dr3(const Vector3d &theta, const Vector3d &x) {
        Matrix3d ddwdr_dr3;
        for (int i = 0; i < 3; ++i)
            ddwdr_dr3.col(i) = get_ddR_i_j(theta, i, 2) * x;
        return ddwdr_dr3;
    }

private:
    constexpr static const double tol = 1e-8;           ///< Tolerance for avoiding singularities of Rodrigues' formula.
    constexpr static const double _PI = 3.14159265359;  ///< Helper.

    /**
     * @brief Get the skew symmetric matrix for a given vector.
     * 
     * @param[in] v The vector.
     * @return The skew symmetric matrix.
     */
    static Matrix3d getSkewSymmetricMatrix(const Vector3d &v) {
        Matrix3d result;
        result << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
        return result;
    }

    /**
     * @brief Safe arcus cosinus.
     * @param[in] x The value.
     * @return safe arcus cosinus.
     */
    static double safeACOS(double x) {
        if (x < -1)
            return _PI;
        if (x > 1)
            return 0;
        return acos(x);
    }

    /**
     * @brief Safe derivative of arcus cosinus.
     * @param[in] x The value.
     * @return Safe derivative of the arcus cosinus.
     */
    static double safeDACOS(double x) {
        if (x < -1 || x > 1)
            return 2 * _PI;
        return -1.0 / sqrt(1.0 - x * x);
    }
};
}  // namespace DCA
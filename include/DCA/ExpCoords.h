#pragma once

#include <DCA/Logger.h>
#include <DCA/Utils.h>
#include <DCA/ddR.h>

#include <array>

namespace DCA {

class ExpCoords {
public:
    //returns exponential map parameterization from rotation matrix
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

    //returns derivative of exponential coorinates with respect to a column of R
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
            dTdRi.col(0) = -0.5 * dacos * theta / theta_norm;
            dTdRi.col(1) = -acos * (theta * theta[2] / theta_norm_3 - Vector3d::UnitZ() / theta_norm);
            dTdRi.col(2) = acos * (theta * theta[1] / theta_norm_3 - Vector3d::UnitY() / theta_norm);
        } else if (index == 1) {
            dTdRi.col(0) = acos * (theta * theta[2] / theta_norm_3 - Vector3d::UnitZ() / theta_norm);
            dTdRi.col(1) = -0.5 * dacos * theta / theta_norm;
            dTdRi.col(2) = -acos * (theta * theta[0] / theta_norm_3 - Vector3d::UnitX() / theta_norm);
        } else if (index == 2) {
            dTdRi.col(0) = -acos * (theta * theta[1] / theta_norm_3 - Vector3d::UnitY() / theta_norm);
            dTdRi.col(1) = acos * (theta * theta[0] / theta_norm_3 - Vector3d::UnitX() / theta_norm);
            dTdRi.col(2) = -0.5 * dacos * theta / theta_norm;
        }
        return dTdRi;
    }

    //returns rotation matrix from exponential map
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

    //returns a Jacobian that tells us how the rotation matrix changes with respect to the exponential map parameters that define it
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

    //returns the derivative of the Jacobian dRi wrt the exponential map parameters at index j
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

    //returns the world coordinates for the vector x that is expressed in local coordinates
    static Vector3d get_w(const Vector3d &theta, const Vector3d &x) {
        Vector3d res = get_R(theta) * x;
        return res;
    }

    //returns a matrix that tells us how w changes wrt the orientation...
    static Matrix3d get_dwdr(const Vector3d &theta, const Vector3d &x) {
        Matrix3d dw_dr;
        for (int i = 0; i < 3; ++i)
            dw_dr.col(i) = get_dRi(theta, i) * x;
        return dw_dr;
    }

    //returns a matrix that tells us how the Jacobian dwdr changes wrt the first component of the orientation...
    static Matrix3d get_ddwdr_dr1(const Vector3d &theta, const Vector3d &x) {
        Matrix3d ddwdr_dr1;
        for (int i = 0; i < 3; ++i)
            ddwdr_dr1.col(i) = get_ddR_i_j(theta, i, 0) * x;
        return ddwdr_dr1;
    }

    //returns a matrix that tells us how the Jacobian dwdr changes wrt the second component of the orientation...
    static Matrix3d get_ddwdr_dr2(const Vector3d &theta, const Vector3d &x) {
        Matrix3d ddwdr_dr2;
        for (int i = 0; i < 3; ++i)
            ddwdr_dr2.col(i) = get_ddR_i_j(theta, i, 1) * x;
        return ddwdr_dr2;
    }

    //returns a matrix that tells us how the Jacobian dwdr changes wrt the third component of the orientation...
    static Matrix3d get_ddwdr_dr3(const Vector3d &theta, const Vector3d &x) {
        Matrix3d ddwdr_dr3;
        for (int i = 0; i < 3; ++i)
            ddwdr_dr3.col(i) = get_ddR_i_j(theta, i, 2) * x;
        return ddwdr_dr3;
    }

private:
    // tolerance for avoiding singularity of Rodrigues' formula
    constexpr static const double tol = 1e-8;
    constexpr static const double PI = 3.14159265359;

    //make skew symmetric for v
    static Matrix3d getSkewSymmetricMatrix(const Vector3d &v) {
        Matrix3d result;
        result << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
        return result;
    }

    static double safeACOS(double x) {
        if (x < -1)
            return PI;
        if (x > 1)
            return 0;
        return acos(x);
    }

    static double safeDACOS(double x) {
        if (x < -1 || x > 1)
            return 2 * PI;
        return 1.0 / sqrt(1.0 - x * x);
    }
};
}  // namespace DCA
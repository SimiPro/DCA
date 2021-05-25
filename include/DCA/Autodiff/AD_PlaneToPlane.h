#pragma once
#include <Eigen/Core>

namespace PlaneToPlaneDistance_CodeGen {
void AD_PlaneToPlaneDistance(const Eigen::Matrix<double, 12, 1>& P,
                             double& objVal);
void AD_PlaneToPlaneDistanceGradient(const Eigen::Matrix<double, 12, 1>& P,
                                     Eigen::Matrix<double, 12, 1>& gradient);
void AD_PlaneToPlaneDistanceHessian(const Eigen::Matrix<double, 12, 1>& P,
                                    Eigen::Matrix<double, 12, 12>& hessian);
}  // namespace PlaneToPlaneDistance_CodeGen

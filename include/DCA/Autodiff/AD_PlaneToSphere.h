#pragma once
#include <Eigen/Core>

namespace PlaneToSphereDistance_CodeGen {
void AD_PlaneToSphereDistance(const Eigen::Matrix<double, 9, 1>& P,
                              const Eigen::Matrix<double, 1, 1>& props,
                              double& objVal);
void AD_PlaneToSphereDistanceGradient(const Eigen::Matrix<double, 9, 1>& P,
                                      const Eigen::Matrix<double, 1, 1>& props,
                                      Eigen::Matrix<double, 9, 1>& gradient);
void AD_PlaneToSphereDistanceHessian(const Eigen::Matrix<double, 9, 1>& P,
                                     const Eigen::Matrix<double, 1, 1>& props,
                                     Eigen::Matrix<double, 9, 9>& hessian);
}  // namespace PlaneToSphereDistance_CodeGen

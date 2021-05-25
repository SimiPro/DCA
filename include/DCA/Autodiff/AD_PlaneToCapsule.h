#pragma once

#include <Eigen/Core>

namespace PlaneToCapsuleDistance_CodeGen {
void AD_PlaneToCapsuleDistance(const Eigen::Matrix<double, 12, 1>& P,
                               const Eigen::Matrix<double, 1, 1>& props,
                               double curveScale, double sigScale,
                               double& objVal);

void AD_PlaneToCapsuleDistanceGradient(const Eigen::Matrix<double, 12, 1>& P,
                                       const Eigen::Matrix<double, 1, 1>& props,
                                       double curveScale, double sigScale,
                                       Eigen::Matrix<double, 12, 1>& gradient);
void AD_PlaneToCapsuleDistanceHessian(
    const Eigen::Matrix<double, 12, 1>& P,
    const Eigen::Matrix<double, 1, 1>& props, double curveScale,
    double sigScale, Eigen::Matrix<double, 12, 12>& hessian);

}  // namespace PlaneToCapsuleDistance_CodeGen

#include <DCA/Interactions/DCAObjective.h>

#include <iostream>

namespace DCA {
namespace Interactions {

std::pair<Vector3d, Vector3d> DCAObjective::compute_Ps(const VectorXd& s, const VectorXd& t) const {
    return {std::visit([&](const auto& primitive) { return primitive.compute_P(get_s(s, A), get_t(t, A)); }, primitive_A),
            std::visit([&](const auto& primitive) { return primitive.compute_P(get_s(s, B), get_t(t, B)); }, primitive_B)};
}

double DCAObjective::compute_D(const VectorXd& s, const VectorXd& t) const {
    check_inputs(s, t);
    auto [P_A, P_B] = compute_Ps(s, t);
    return (P_A - P_B).squaredNorm();
}

void DCAObjective::compute_pDpT(VectorXd& pDpT, const VectorXd& s, const VectorXd& t) const {
    check_inputs(s, t);
    auto [P_A, P_B] = compute_Ps(s, t);
    auto [SIZE_T_A, SIZE_T_B] = get_sizes_t();

    Eigen::Matrix<double, 3, -1> dPdT_A = std::visit([&](const auto& primitive) { return primitive.compute_dPdT(get_s(s, A), get_t(t, A)); }, primitive_A);
    Eigen::Matrix<double, 3, -1> dPdT_B = std::visit([&](const auto& primitive) { return primitive.compute_dPdT(get_s(s, B), get_t(t, B)); }, primitive_B);

    pDpT.resize(t.size());
    pDpT.segment(0, SIZE_T_A) = 2.0 * dPdT_A.transpose() * (P_A - P_B);
    pDpT.segment(SIZE_T_A, SIZE_T_B) = -2.0 * dPdT_B.transpose() * (P_A - P_B);
}

void DCAObjective::compute_p2DpT2(MatrixXd& p2DpT2, const VectorXd& s, const VectorXd& t) const {
    check_inputs(s, t);
    auto [SIZE_T_A, SIZE_T_B] = get_sizes_t();

    Eigen::Matrix<double, 3, -1> dPdT_A = std::visit([&](const auto& primitive) { return primitive.compute_dPdT(get_s(s, A), get_t(t, A)); }, primitive_A);
    Eigen::Matrix<double, 3, -1> dPdT_B = std::visit([&](const auto& primitive) { return primitive.compute_dPdT(get_s(s, B), get_t(t, B)); }, primitive_B);

    p2DpT2.resize(t.size(), t.size());
    p2DpT2.block(0, 0, SIZE_T_A, SIZE_T_A) = 2.0 * dPdT_A.transpose() * dPdT_A;
    p2DpT2.block(SIZE_T_A, 0, SIZE_T_B, SIZE_T_A) = -2.0 * dPdT_B.transpose() * dPdT_A;
    p2DpT2.block(0, SIZE_T_A, SIZE_T_A, SIZE_T_B) = -2.0 * dPdT_A.transpose() * dPdT_B;
    p2DpT2.block(SIZE_T_A, SIZE_T_A, SIZE_T_B, SIZE_T_B) = 2.0 * dPdT_B.transpose() * dPdT_B;
}

void DCAObjective::compute_pDpS(VectorXd& pDpS, const VectorXd& s, const VectorXd& t) const {
    check_inputs(s, t);
    auto [P_A, P_B] = compute_Ps(s, t);

    Eigen::Matrix<double, 3, 6> dPdS_A = std::visit([&](const auto& primitive) { return primitive.compute_dPdS(get_s(s, A), get_t(t, A)); }, primitive_A);
    Eigen::Matrix<double, 3, 6> dPdS_B = std::visit([&](const auto& primitive) { return primitive.compute_dPdS(get_s(s, B), get_t(t, B)); }, primitive_B);

    pDpS.resize(s.size());
    pDpS.segment<6>(0) = 2.0 * dPdS_A.transpose() * (P_A - P_B);
    pDpS.segment<6>(6) = -2.0 * dPdS_B.transpose() * (P_A - P_B);
}

void DCAObjective::compute_p2DpS2(MatrixXd& p2DpS2, const VectorXd& s, const VectorXd& t) const {
    check_inputs(s, t);

    auto [P_A, P_B] = compute_Ps(s, t);

    Eigen::Matrix<double, 3, 6> dPdS_A = std::visit([&](const auto& primitive) { return primitive.compute_dPdS(get_s(s, A), get_t(t, A)); }, primitive_A);
    Eigen::Matrix<double, 3, 6> dPdS_B = std::visit([&](const auto& primitive) { return primitive.compute_dPdS(get_s(s, B), get_t(t, B)); }, primitive_B);

    std::array<Eigen::Matrix<double, 3, 6>, 6> d2PdS2_A =
        std::visit([&](const auto& primitive) { return primitive.compute_d2PdS2(get_s(s, A), get_t(t, A)); }, primitive_A);
    std::array<Eigen::Matrix<double, 3, 6>, 6> d2PdS2_B =
        std::visit([&](const auto& primitive) { return primitive.compute_d2PdS2(get_s(s, B), get_t(t, B)); }, primitive_B);

    p2DpS2.resize(s.size(), s.size());

    p2DpS2.block<6, 6>(0, 0) = 2.0 * dPdS_A.transpose() * dPdS_A;
    p2DpS2.block<6, 6>(0, 6) = -2.0 * dPdS_A.transpose() * dPdS_B;
    p2DpS2.block<6, 6>(6, 0) = -2.0 * dPdS_B.transpose() * dPdS_A;
    p2DpS2.block<6, 6>(6, 6) = 2.0 * dPdS_B.transpose() * dPdS_B;

    for (int i = 0; i < 6; i++) {
        p2DpS2.block<6, 1>(0, i) += 2.0 * d2PdS2_A[i].transpose() * (P_A - P_B);
        p2DpS2.block<6, 1>(6, 6 + i) -= 2.0 * d2PdS2_B[i].transpose() * (P_A - P_B);
    }
}

void DCAObjective::compute_p2DpTpS(MatrixXd& p2DpTpS, const VectorXd& s, const VectorXd& t) const {
    check_inputs(s, t);
    auto [P_A, P_B] = compute_Ps(s, t);
    auto [SIZE_T_A, SIZE_T_B] = get_sizes_t();

    Eigen::Matrix<double, 3, -1> dPdT_A = std::visit([&](const auto& primitive) { return primitive.compute_dPdT(get_s(s, A), get_t(t, A)); }, primitive_A);
    Eigen::Matrix<double, 3, -1> dPdT_B = std::visit([&](const auto& primitive) { return primitive.compute_dPdT(get_s(s, B), get_t(t, B)); }, primitive_B);

    Eigen::Matrix<double, 3, 6> dPdS_A = std::visit([&](const auto& primitive) { return primitive.compute_dPdS(get_s(s, A), get_t(t, A)); }, primitive_A);
    Eigen::Matrix<double, 3, 6> dPdS_B = std::visit([&](const auto& primitive) { return primitive.compute_dPdS(get_s(s, B), get_t(t, B)); }, primitive_B);

    std::vector<Eigen::Matrix<double, 3, 6>> d2PdSdT_A =
        std::visit([&](const auto& primitive) { return primitive.compute_d2PdSdT(get_s(s, A), get_t(t, A)); }, primitive_A);
    std::vector<Eigen::Matrix<double, 3, 6>> d2PdSdT_B =
        std::visit([&](const auto& primitive) { return primitive.compute_d2PdSdT(get_s(s, B), get_t(t, B)); }, primitive_B);

    p2DpTpS.resize(t.size(), s.size());

    for (int a = 0; a < SIZE_T_A; a++) {
        Vector3d dPdT_A_part = dPdT_A.block(0, a, 3, 1);
        p2DpTpS.block(a, 0, 1, 6) = 2.0 * ((P_A - P_B).transpose() * d2PdSdT_A[a] + dPdT_A_part.transpose() * dPdS_A);
        p2DpTpS.block(a, 6, 1, 6) = -2.0 * dPdT_A_part.transpose() * dPdS_B;
    }

    for (int b = 0; b < SIZE_T_B; b++) {
        Vector3d dPdT_B_part = dPdT_B.block(0, b, 3, 1);
        p2DpTpS.block(SIZE_T_A + b, 0, 1, 6) = -2.0 * dPdT_B_part.transpose() * dPdS_A;
        p2DpTpS.block(SIZE_T_A + b, 6, 1, 6) = 2.0 * (-(P_A - P_B).transpose() * d2PdSdT_B[b] + dPdT_B_part.transpose() * dPdS_B);
    }
}

double DCAObjective::compute_O(const VectorXd& s, const VectorXd& t) const {
    double value = 0.0;

    //--- Shortest distance
    value += compute_D(s, t);

    //--- Regularizer
    value += regularizerWeight * 0.5 * t.dot(t);

    //--- Constraints
    for (int i = 0; i < (int)t.size(); i++) {
        value += constraintWeight * unilateralConstraint.evaluate(-t[i] - 1.0);
        value += constraintWeight * unilateralConstraint.evaluate(t[i] - 1.0);
    }

    return value;
}

void DCAObjective::add_pOpT_To(VectorXd& pOpT, const VectorXd& s, const VectorXd& t) const {
    //--- Shortest distance
    VectorXd pDpT(t.size());
    compute_pDpT(pDpT, s, t);
    pOpT += pDpT;

    //--- Regularizer
    pOpT += regularizerWeight * t;

    //--- Constraint
    for (int i = 0; i < (int)t.size(); i++) {
        pOpT[i] -= constraintWeight * unilateralConstraint.computeDerivative(-t[i] - 1.0);
        pOpT[i] += constraintWeight * unilateralConstraint.computeDerivative(t[i] - 1.0);
    }
}

void DCAObjective::add_pOpS_To(VectorXd& pOpS, const VectorXd& s, const VectorXd& t) const {
    //--- Shortest distance
    VectorXd pDpS(s.size());
    compute_pDpS(pDpS, s, t);
    pOpS += pDpS;
}

void DCAObjective::add_p2OpT2_To(MatrixXd& p2OpT2, const VectorXd& s, const VectorXd& t) const {
    //--- Shortest distance
    MatrixXd p2DpT2(t.size(), t.size());
    compute_p2DpT2(p2DpT2, s, t);
    p2OpT2 += p2DpT2;

    //--- Regularizer
    VectorXd reg(t.size());
    reg.setConstant(regularizerWeight);
    p2OpT2 += reg.asDiagonal();

    //--- Constraint
    for (int i = 0; i < (int)t.size(); i++) {
        p2OpT2(i, i) += constraintWeight * unilateralConstraint.computeSecondDerivative(-t[i] - 1.0);
        p2OpT2(i, i) += constraintWeight * unilateralConstraint.computeSecondDerivative(t[i] - 1.0);
    }
}

void DCAObjective::add_p2OpS2_To(MatrixXd& p2OpS2, const VectorXd& s, const VectorXd& t) const {
    //--- Shortest distance
    MatrixXd p2DpS2(s.size(), s.size());
    compute_p2DpS2(p2DpS2, s, t);
    p2OpS2 += p2DpS2;
}

void DCAObjective::add_p2OpTpS_To(MatrixXd& p2OpTpS, const VectorXd& s, const VectorXd& t) const {
    //--- Shortest distance
    MatrixXd p2DpTpS(t.size(), s.size());
    compute_p2DpTpS(p2DpTpS, s, t);
    p2OpTpS += p2DpTpS;
}

std::pair<int, int> DCAObjective::get_sizes_t() const {
    return {std::visit([&](const auto& primitive) { return primitive.SIZE_T(); }, primitive_A),
            std::visit([&](const auto& primitive) { return primitive.SIZE_T(); }, primitive_B)};
}

int DCAObjective::get_total_size_t() const {
    auto [SIZE_T_A, SIZE_T_B] = get_sizes_t();
    return SIZE_T_A + SIZE_T_B;
}

Vector6d DCAObjective::get_s(const VectorXd& s, PRIMITIVE P) const {
    if (P == A)
        return s.segment(0, 6);
    else if (P == B)
        return s.segment(6, 6);
    throw std::runtime_error("ERROR in DCAObjective::get_s -> something is wrong with s...");
    return Vector6d::Zero();
}

VectorXd DCAObjective::get_t(const VectorXd& t, PRIMITIVE P) const {
    auto [SIZE_T_A, SIZE_T_B] = get_sizes_t();
    if (P == A)
        return t.segment(0, SIZE_T_A);
    else if (P == B)
        return t.segment(SIZE_T_A, SIZE_T_B);
    throw std::runtime_error("ERROR in DCAObjective::get_t -> something is wrong with t...");
    return VectorXd::Zero(0);
}

void DCAObjective::check_inputs(const VectorXd& s, const VectorXd& t) const {
    if ((int)s.size() != 12)
        throw std::runtime_error("ERROR in DCAObjective::check_inputs -> invalid input s!");

    if ((int)t.size() != get_total_size_t())
        throw std::runtime_error("ERROR in DCAObjective::check_inputs -> invalid input t!");
}

std::string DCAObjective::get_primitives_description(primitive_t primitive_A, primitive_t primitive_B) {
    return std::visit([&](const auto& primitive) { return primitive.description; }, primitive_A) + "VS" +
           std::visit([&](const auto& primitive) { return primitive.description; }, primitive_B);
}
}  // namespace Interactions
}  // namespace DCA
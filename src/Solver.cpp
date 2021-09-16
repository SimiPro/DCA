#include <DCA/Interactions/Solver.h>
#include <DCA/Utils/Logger.h>

namespace DCA {
namespace Interactions {

Solver::Solver(primitive_t primitive_A, primitive_t primitive_B)
    : Opt::FiniteDifference(DCAObjective::get_primitives_description(primitive_A, primitive_B) + "-Solver"), objective(primitive_A, primitive_B) {
    deltaFD = 1e-8;
    absTol = 1e-4;
}

void Solver::compute_t(VectorXd& t, const VectorXd& s, bool forFD) const {
    if ((int)t.size() != objective.get_total_size_t())
        throw std::runtime_error("Solver::compute_t: input size of t is not correct!");

    optimizer.solverResidual = forFD ? 1e-12 : 1e-6;
    optimizer.optimize(t, s, objective, 1000);
}

void Solver::compute_dtds(MatrixXd& dtds, const VectorXd& s, const VectorXd& t) const {
    MatrixXd p2OpT2(t.size(), t.size()), p2OpTpS(t.size(), s.size());
    objective.get_p2OpT2(p2OpT2, s, t);
    objective.get_p2OpTpS(p2OpTpS, s, t);

    //Solve linear system
    Eigen::LLT<Eigen::MatrixXd, Eigen::Lower> solver;
    solver.compute(p2OpT2.triangularView<Eigen::Lower>());
    dtds = solver.solve(-p2OpTpS);
    if (solver.info() != Eigen::Success)
        Logger::print(Logger::RED, "WARNING -> Solver::compute_dtds: Y = A*X (dense) solve unsuccessful!\n");
}

void Solver::test_dtds_WithFD(const VectorXd& s, const VectorXd& t) const {
    auto eval = [&](VectorXd& vec, const VectorXd& s, const VectorXd& t) -> void {
        vec = t;
        compute_t(vec, s, true);
    };
    auto anal = [&](MatrixXd& mat, const VectorXd& s, const VectorXd& t) -> void { compute_dtds(mat, s, t); };
    testMatrix(eval, anal, s, t, D_S, "dtds", (int)t.size());
}

double Solver::compute_D(const VectorXd& s, const VectorXd& t) const {
    are_t_and_s_in_sync(s, t);
    double r_A = std::visit([&](const auto& primitive) { return primitive.radius; }, objective.primitive_A);
    double r_B = std::visit([&](const auto& primitive) { return primitive.radius; }, objective.primitive_B);
    return objective.compute_D(s, t) - (r_A + r_B) * (r_A + r_B);
}

void Solver::compute_dDdS(VectorXd& dDdS, const VectorXd& s, const VectorXd& t) const {
    are_t_and_s_in_sync(s, t);

    VectorXd pDpT(t.size()), pDpS(s.size());
    objective.compute_pDpT(pDpT, s, t);
    objective.compute_pDpS(pDpS, s, t);

    MatrixXd dtds(t.size(), s.size());
    compute_dtds(dtds, s, t);

    dDdS = dtds.transpose() * pDpT + pDpS;
}

void Solver::test_dDdS_WithFD(const VectorXd& s, const VectorXd& t) const {
    auto eval = [&](const VectorXd& s, const VectorXd& t) -> double {
        VectorXd t_tmp(t);
        compute_t(t_tmp, s, true);
        return compute_D(s, t_tmp);
    };
    auto anal = [&](VectorXd& vec, const VectorXd& s, const VectorXd& t) -> void { compute_dDdS(vec, s, t); };
    testVector(eval, anal, s, t, D_S, "dDdS");
}

void Solver::compute_d2DdS2(MatrixXd& d2DdS2, const VectorXd& s, const VectorXd& t) const {
    are_t_and_s_in_sync(s, t);

    MatrixXd p2DpT2(t.size(), t.size()), p2DpS2(s.size(), s.size()), p2DpTpS(t.size(), s.size());
    objective.compute_p2DpT2(p2DpT2, s, t);
    objective.compute_p2DpS2(p2DpS2, s, t);
    objective.compute_p2DpTpS(p2DpTpS, s, t);

    MatrixXd dtds(t.size(), s.size());
    compute_dtds(dtds, s, t);

    d2DdS2 = (dtds.transpose() * p2DpT2 + 2.0 * p2DpTpS.transpose()) * dtds + p2DpS2;
}

void Solver::test_d2DdS2_WithFD(const VectorXd& s, const VectorXd& t) const {
    auto eval = [&](VectorXd& vec, const VectorXd& s, const VectorXd& t) -> void {
        VectorXd t_tmp(t);
        compute_t(t_tmp, s, true);
        compute_dDdS(vec, s, t_tmp);
    };
    auto anal = [&](MatrixXd& mat, const VectorXd& s, const VectorXd& t) -> void { compute_d2DdS2(mat, s, t); };
    testMatrix(eval, anal, s, t, D_S, "d2DdS2", (int)s.size());
}

bool Solver::are_t_and_s_in_sync(const VectorXd& s, const VectorXd& t) const {
    if ((int)t.size() != objective.get_total_size_t())
        throw std::runtime_error("Solver::are_t_and_s_in_sync: input size of t is not correct!");

    VectorXd grad(t.size());
    objective.get_pOpT(grad, s, t);
    if (grad.norm() > 1e-6) {
        Logger::print(Logger::RED, "ERROR in Solver::are_t_and_s_in_sync -> s and t are NOT in sync! Gradient norm is: %lf", grad.norm());
        throw std::runtime_error("Solver::are_t_and_s_in_sync -> s and t are NOT in sync!");
        return false;
    }
    return true;
}

std::pair<Vector3d, Vector3d> Solver::compute_closest_points(const VectorXd& s, const VectorXd& t) const {
    are_t_and_s_in_sync(s, t);
    return objective.compute_Ps(s, t);
}

Vector12d Solver::get_s_from_primitives() const {
    Vector6d s_A = std::visit([&](const auto& primitive) { return primitive.get_s(); }, objective.primitive_A);
    Vector6d s_B = std::visit([&](const auto& primitive) { return primitive.get_s(); }, objective.primitive_B);
    Vector12d s;
    s << s_A, s_B;
    return s;
}

int Solver::get_total_size_t() const {
    return objective.get_total_size_t();
}

void Solver::test_derivatives() const {
    VectorXd s = get_s_from_primitives();
    VectorXd t = VectorXd::Zero(objective.get_total_size_t());
    compute_t(t, s, true);

    objective.test_pOpT_WithFD(s, t);
    objective.test_p2OpT2_WithFD(s, t);

    objective.test_pOpS_WithFD(s, t);
    objective.test_p2OpS2_WithFD(s, t);
    objective.test_p2OpTpS_WithFD(s, t);

    test_dtds_WithFD(s, t);
    test_dDdS_WithFD(s, t);
    //test_d2DdS2_WithFD(s, t);
}

}  // namespace Interactions
}  // namespace DCA
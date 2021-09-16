#include <DCA/API.h>
#include <DCA/Interactions/Solver.h>

namespace DCA {

double API::compute_D(const primitive_t& p_a, const primitive_t& p_b) {
    Interactions::Solver solver(p_a, p_b);
    VectorXd s = solver.get_s_from_primitives();
    VectorXd t = VectorXd::Zero(solver.get_total_size_t());
    solver.compute_t(t, s, false);
    return solver.compute_D(s, t);
}

void API::compute_dDdS(Vector12d& dDdS, const primitive_t& p_a, const primitive_t& p_b) {
    Interactions::Solver solver(p_a, p_b);
    VectorXd s = solver.get_s_from_primitives();
    VectorXd t = VectorXd::Zero(solver.get_total_size_t());
    solver.compute_t(t, s, false);
    VectorXd dDdS_tmp(12);
    solver.compute_dDdS(dDdS_tmp, s, t);
    dDdS = dDdS_tmp;
}

void API::compute_d2DdS2(Matrix12d& d2DdS2, const primitive_t& p_a, const primitive_t& p_b) {
    Interactions::Solver solver(p_a, p_b);
    VectorXd s = solver.get_s_from_primitives();
    VectorXd t = VectorXd::Zero(solver.get_total_size_t());
    solver.compute_t(t, s, false);
    MatrixXd d2DdS2_tmp(12, 12);
    solver.compute_d2DdS2(d2DdS2_tmp, s, t);
    d2DdS2 = d2DdS2_tmp;
}

double API::compute_D(const primitive_t& p_a, const primitive_t& p_b, const VectorXd& t) {
    Interactions::Solver solver(p_a, p_b);
    return solver.compute_D(solver.get_s_from_primitives(), t);
}

void API::compute_dDdS(Vector12d& dDdS, const primitive_t& p_a, const primitive_t& p_b, const VectorXd& t) {
    Interactions::Solver solver(p_a, p_b);
    VectorXd dDdS_tmp(12);
    solver.compute_dDdS(dDdS_tmp, solver.get_s_from_primitives(), t);
    dDdS = dDdS_tmp;
}

void API::compute_d2DdS2(Matrix12d& d2DdS2, const primitive_t& p_a, const primitive_t& p_b, const VectorXd& t) {
    Interactions::Solver solver(p_a, p_b);
    MatrixXd d2DdS2_tmp(12, 12);
    solver.compute_d2DdS2(d2DdS2_tmp, solver.get_s_from_primitives(), t);
    d2DdS2 = d2DdS2_tmp;
}

std::pair<Vector3d, Vector3d> API::compute_closest_points(const primitive_t& p_a, const primitive_t& p_b) {
    Interactions::Solver solver(p_a, p_b);
    VectorXd s = solver.get_s_from_primitives();
    VectorXd t = VectorXd::Zero(solver.get_total_size_t());
    solver.compute_t(t, s, false);
    return solver.compute_closest_points(s, t);
}

std::pair<Vector3d, Vector3d> API::compute_closest_points(const primitive_t& p_a, const primitive_t& p_b, const VectorXd& t) {
    Interactions::Solver solver(p_a, p_b);
    return solver.compute_closest_points(solver.get_s_from_primitives(), t);
}

void API::compute_t(VectorXd& t, const primitive_t& p_a, const primitive_t& p_b) {
    Interactions::Solver solver(p_a, p_b);
    t.resize(solver.get_total_size_t());
    t.setZero();
    solver.compute_t(t, solver.get_s_from_primitives(), false);
}
}  // namespace DCA
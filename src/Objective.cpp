#include <DCA/Opt/Objective.h>

namespace DCA {
namespace Opt {

void Objective::add_pOpT_To(VectorXd& pOpT, const VectorXd& s, const VectorXd& t) const {
    auto eval = [&](const VectorXd& s, const VectorXd& t) -> double { return compute_O(s, t); };
    VectorXd vec(t.size());
    estimateVector(vec, eval, s, t, D_T);
    pOpT += vec;
}

void Objective::add_pOpS_To(VectorXd& pOpS, const VectorXd& s, const VectorXd& t) const {
    auto eval = [&](const VectorXd& s, const VectorXd& t) -> double { return compute_O(s, t); };
    VectorXd vec(s.size());
    estimateVector(vec, eval, s, t, D_S);
    pOpS += vec;
}

void Objective::add_p2OpT2_To(MatrixXd& p2OpT2, const VectorXd& s, const VectorXd& t) const {
    auto eval = [&](VectorXd& vec, const VectorXd& s, const VectorXd& t) -> void { get_pOpT(vec, s, t); };
    MatrixXd mat(t.size(), t.size());
    estimateMatrix(mat, eval, s, t, D_T, (int)t.size());
    p2OpT2 += mat;
}

void Objective::add_p2OpS2_To(MatrixXd& p2OpS2, const VectorXd& s, const VectorXd& t) const {
    auto eval = [&](VectorXd& vec, const VectorXd& s, const VectorXd& t) -> void { get_pOpS(vec, s, t); };
    MatrixXd mat(s.size(), s.size());
    estimateMatrix(mat, eval, s, t, D_S, (int)s.size());
    p2OpS2 += mat;
}

void Objective::add_p2OpTpS_To(MatrixXd& p2OpTpS, const VectorXd& s, const VectorXd& t) const {
    auto eval = [&](VectorXd& vec, const VectorXd& s, const VectorXd& t) -> void { get_pOpT(vec, s, t); };
    MatrixXd mat(t.size(), s.size());
    estimateMatrix(mat, eval, s, t, D_S, (int)t.size());
    p2OpTpS += mat;
}

void Objective::get_pOpT(VectorXd& pOpT, const VectorXd& s, const VectorXd& t) const {
    pOpT.resize(t.size());
    pOpT.setZero();
    add_pOpT_To(pOpT, s, t);
}

void Objective::get_pOpS(VectorXd& pOpS, const VectorXd& s, const VectorXd& t) const {
    pOpS.resize(s.size());
    pOpS.setZero();
    add_pOpS_To(pOpS, s, t);
}

void Objective::get_p2OpT2(MatrixXd& p2OpT2, const VectorXd& s, const VectorXd& t) const {
    p2OpT2.resize(t.size(), t.size());
    p2OpT2.setZero();
    add_p2OpT2_To(p2OpT2, s, t);
}

void Objective::get_p2OpS2(MatrixXd& p2OpS2, const VectorXd& s, const VectorXd& t) const {
    p2OpS2.resize(s.size(), s.size());
    p2OpS2.setZero();
    add_p2OpS2_To(p2OpS2, s, t);
}

void Objective::get_p2OpTpS(MatrixXd& p2OpTpS, const VectorXd& s, const VectorXd& t) const {
    p2OpTpS.resize(t.size(), s.size());
    p2OpTpS.setZero();
    add_p2OpTpS_To(p2OpTpS, s, t);
}

void Objective::test_pOpT_WithFD(const VectorXd& s, const VectorXd& t) const {
    auto eval = [&](const VectorXd& s, const VectorXd& t) -> double { return compute_O(s, t); };
    auto anal = [&](VectorXd& vec, const VectorXd& s, const VectorXd& t) -> void { get_pOpT(vec, s, t); };
    testVector(eval, anal, s, t, D_T, "pOpT");
}

void Objective::test_pOpS_WithFD(const VectorXd& s, const VectorXd& t) const {
    auto eval = [&](const VectorXd& s, const VectorXd& t) -> double { return compute_O(s, t); };
    auto anal = [&](VectorXd& vec, const VectorXd& s, const VectorXd& t) -> void { get_pOpS(vec, s, t); };
    testVector(eval, anal, s, t, D_S, "pOpS");
}

void Objective::test_p2OpT2_WithFD(const VectorXd& s, const VectorXd& t) const {
    auto eval = [&](VectorXd& vec, const VectorXd& s, const VectorXd& t) -> void { get_pOpT(vec, s, t); };
    auto anal = [&](MatrixXd& mat, const VectorXd& s, const VectorXd& t) -> void { get_p2OpT2(mat, s, t); };
    testMatrix(eval, anal, s, t, D_T, "p2OpT2", (int)t.size());
}

void Objective::test_p2OpS2_WithFD(const VectorXd& s, const VectorXd& t) const {
    auto eval = [&](VectorXd& vec, const VectorXd& s, const VectorXd& t) -> void { get_pOpS(vec, s, t); };
    auto anal = [&](MatrixXd& mat, const VectorXd& s, const VectorXd& t) -> void { get_p2OpS2(mat, s, t); };
    testMatrix(eval, anal, s, t, D_S, "p2OpS2", (int)s.size());
}

void Objective::test_p2OpTpS_WithFD(const VectorXd& s, const VectorXd& t) const {
    auto eval = [&](VectorXd& vec, const VectorXd& s, const VectorXd& t) -> void { get_pOpT(vec, s, t); };
    auto anal = [&](MatrixXd& mat, const VectorXd& s, const VectorXd& t) -> void { get_p2OpTpS(mat, s, t); };
    testMatrix(eval, anal, s, t, D_S, "p2OpTpS", (int)t.size());
}

}  // namespace Opt
}  // namespace DCA
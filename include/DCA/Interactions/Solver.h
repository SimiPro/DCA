#pragma once

#include <DCA/Interactions/DCAObjective.h>
#include <DCA/Opt/NewtonOptimizer.h>

namespace DCA {
namespace Interactions {

class Solver : public Opt::FiniteDifference {
public:
    Solver(primitive_t primitive_A, primitive_t primitive_B);
    ~Solver() {}

    void compute_t(VectorXd& t, const VectorXd& s, bool forFD = false) const;

    void compute_dtds(MatrixXd& dtds, const VectorXd& s, const VectorXd& t) const;
    void test_dtds_WithFD(const VectorXd& s, const VectorXd& t) const;

    double compute_D(const VectorXd& s, const VectorXd& t) const;

    void compute_dDdS(VectorXd& dDdS, const VectorXd& s, const VectorXd& t) const;
    void test_dDdS_WithFD(const VectorXd& s, const VectorXd& t) const;

    void compute_d2DdS2(MatrixXd& d2DdS2, const VectorXd& s, const VectorXd& t) const;
    void test_d2DdS2_WithFD(const VectorXd& s, const VectorXd& t) const;

    std::pair<Vector3d, Vector3d> compute_closest_points(const VectorXd& s, const VectorXd& t) const;

    bool are_t_and_s_in_sync(const VectorXd& s, const VectorXd& t) const;
    Vector12d get_s_from_primitives() const;
    int get_total_size_t() const;
    void test_derivatives() const;

private:
    DCAObjective objective;
    mutable Opt::NewtonOptimizer optimizer;
};

}  // namespace Interactions

}  // namespace DCA
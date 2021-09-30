#pragma once

#include <DCA/FiniteDifference.h>
#include <DCA/Primitives.h>
#include <DCA/SoftUnilateralConstraint.h>

namespace DCA {

class Objective : public FiniteDifference {
public:
    Objective(primitive_t primitive_A, primitive_t primitive_B)
        : FiniteDifference(get_primitives_description(primitive_A, primitive_B) + "-Objective"), primitive_A(primitive_A), primitive_B(primitive_B) {}
    ~Objective() {}

    //--- Compute points on primitives based on s and t
    std::pair<Vector3d, Vector3d> compute_Ps(const VectorXd& s, const VectorXd& t) const;

    //--- Compute D and its derivatives
    double compute_D(const VectorXd& s, const VectorXd& t) const;
    void compute_pDpT(VectorXd& pDpT, const VectorXd& s, const VectorXd& t) const;
    void compute_p2DpT2(MatrixXd& p2DpT2, const VectorXd& s, const VectorXd& t) const;
    void compute_pDpS(VectorXd& pDpS, const VectorXd& s, const VectorXd& t) const;
    void compute_p2DpS2(MatrixXd& p2DpS2, const VectorXd& s, const VectorXd& t) const;
    void compute_p2DpTpS(MatrixXd& p2DpTpS, const VectorXd& s, const VectorXd& t) const;

    //--- Compute O and its derivatives
    double compute_O(const VectorXd& s, const VectorXd& t) const;
    void compute_pOpT(VectorXd& pOpT, const VectorXd& s, const VectorXd& t) const;
    void compute_pOpS(VectorXd& pOpS, const VectorXd& s, const VectorXd& t) const;
    void compute_p2OpT2(MatrixXd& p2OpT2, const VectorXd& s, const VectorXd& t) const;
    void compute_p2OpS2(MatrixXd& p2OpS2, const VectorXd& s, const VectorXd& t) const;
    void compute_p2OpTpS(MatrixXd& p2OpTpS, const VectorXd& s, const VectorXd& t) const;

    //--- Finite difference tests
    void test_pOpT_WithFD(const VectorXd& s, const VectorXd& t) const;
    void test_pOpS_WithFD(const VectorXd& s, const VectorXd& t) const;
    void test_p2OpT2_WithFD(const VectorXd& s, const VectorXd& t) const;
    void test_p2OpS2_WithFD(const VectorXd& s, const VectorXd& t) const;
    void test_p2OpTpS_WithFD(const VectorXd& s, const VectorXd& t) const;

    //--- Public helpers
    std::pair<int, int> get_sizes_t() const;
    int get_total_size_t() const;
    static std::string get_primitives_description(primitive_t primitive_A, primitive_t primitive_B);

private:
    //--- Private helpers
    enum PRIMITIVE { A, B };
    Vector6d get_s(const VectorXd& s, PRIMITIVE P) const;
    VectorXd get_t(const VectorXd& t, PRIMITIVE P) const;
    void check_inputs(const VectorXd& s, const VectorXd& t) const;

public:
    double regularizerWeight = 0.01;
    double constraintWeight = 10.0;

    SoftUnilateralConstraint unilateralConstraint = SoftUnilateralConstraint(0.0, 1.0, 1e-3);

    const primitive_t primitive_A, primitive_B;
};

}  // namespace DCA
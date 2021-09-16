#pragma once

#include <DCA/Opt/Objective.h>
#include <DCA/Opt/SoftUnilateralConstraint.h>
#include <DCA/Utils/Primitives.h>

namespace DCA {
namespace Interactions {

class DCAObjective : public Opt::Objective {
public:
    DCAObjective(primitive_t primitive_A, primitive_t primitive_B)
        : Opt::Objective(get_primitives_description(primitive_A, primitive_B) + "-Objective"), primitive_A(primitive_A), primitive_B(primitive_B) {}
    ~DCAObjective() {}

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
    double compute_O(const VectorXd& s, const VectorXd& t) const override;
    void add_pOpT_To(VectorXd& pOpT, const VectorXd& s, const VectorXd& t) const override;
    void add_pOpS_To(VectorXd& pOpS, const VectorXd& s, const VectorXd& t) const override;
    void add_p2OpT2_To(MatrixXd& p2OpT2, const VectorXd& s, const VectorXd& t) const override;
    void add_p2OpS2_To(MatrixXd& p2OpS2, const VectorXd& s, const VectorXd& t) const override;
    void add_p2OpTpS_To(MatrixXd& p2OpTpS, const VectorXd& s, const VectorXd& t) const override;

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
    double regularizerWeight = 0.01; // 0.01
    double constraintWeight = 10.0;

    Opt::SoftUnilateralConstraint unilateralConstraint = Opt::SoftUnilateralConstraint(0.0, 1.0, 1e-3);

    const primitive_t primitive_A, primitive_B;
};

}  // namespace Interactions
}  // namespace DCA
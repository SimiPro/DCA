#pragma once

#include <DCA/Opt/FiniteDifference.h>

namespace DCA {
namespace Opt {

class Objective : public FiniteDifference {
public:
    Objective(std::string description) : FiniteDifference(description) {}
    virtual ~Objective() {}

    //--- Objective and derivatives (with FD as default) ---
    virtual double compute_O(const VectorXd& s, const VectorXd& t) const = 0;

    virtual void add_pOpT_To(VectorXd& pOpT, const VectorXd& s, const VectorXd& t) const;
    virtual void add_pOpS_To(VectorXd& pOpS, const VectorXd& s, const VectorXd& t) const;

    virtual void add_p2OpT2_To(MatrixXd& p2OpT2, const VectorXd& s, const VectorXd& t) const;
    virtual void add_p2OpS2_To(MatrixXd& p2OpS2, const VectorXd& s, const VectorXd& t) const;
    virtual void add_p2OpTpS_To(MatrixXd& p2OpTpS, const VectorXd& s, const VectorXd& t) const;

    void get_pOpT(VectorXd& pOpT, const VectorXd& s, const VectorXd& t) const;
    void get_pOpS(VectorXd& pOpS, const VectorXd& s, const VectorXd& t) const;

    void get_p2OpT2(MatrixXd& p2OpT2, const VectorXd& s, const VectorXd& t) const;
    void get_p2OpS2(MatrixXd& p2OpS2, const VectorXd& s, const VectorXd& t) const;
    void get_p2OpTpS(MatrixXd& p2OpTpS, const VectorXd& s, const VectorXd& t) const;

    //--- Test derivatives with finite difference ---
    virtual void test_pOpT_WithFD(const VectorXd& s, const VectorXd& t) const;
    virtual void test_pOpS_WithFD(const VectorXd& s, const VectorXd& t) const;

    virtual void test_p2OpT2_WithFD(const VectorXd& s, const VectorXd& t) const;
    virtual void test_p2OpS2_WithFD(const VectorXd& s, const VectorXd& t) const;
    virtual void test_p2OpTpS_WithFD(const VectorXd& s, const VectorXd& t) const;
};

}  // namespace Opt
}  // namespace DCA
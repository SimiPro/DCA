#pragma once

#include <DCA/FiniteDifference.h>
#include <DCA/Primitives.h>
#include <DCA/SoftUnilateralConstraint.h>

namespace DCA {

/**
 * @brief The main objective which used for computing distances and their derivatives.
 */
class Objective : public FiniteDifference {
public:
    /**
     * @brief Create an objective for a pair of primitives.
     * @param[in] primitive_A The first primitive.
     * @param[in] primitive_B The second primitive.
     */
    Objective(primitive_t primitive_A, primitive_t primitive_B)
        : FiniteDifference(get_primitives_description(primitive_A, primitive_B) + "-Objective"), primitive_A(primitive_A), primitive_B(primitive_B) {}

    /**
     * @brief Default deconstructor.
     */
    ~Objective() = default;

    /**
     * @brief Compute points on the primitives based on s and t.
     * @param[in] s The state of both primitives, stacked.
     * @param[in] t The parameterization of both primitives, stacked.
     */
    std::pair<Vector3d, Vector3d> compute_Ps(const VectorXd& s, const VectorXd& t) const;

    /**
     * @name D
     * @brief D and its derivatives
     */
    /** @{ */  // start of group

    /**
     * @brief Compute the distance.
     * @param[in] s The state of the two primitives.
     * @param[in] t The parameterization of the two primitives.
     * @return The distance between the two primitives.
     */
    double compute_D(const VectorXd& s, const VectorXd& t) const;

    /**
     * @brief Compute the partial derivative of the distance w.r.t. t.
     * @param[out] pDpT \f$ \frac{\partial D}{\partial t} \f$
     * @param[in] s The state of the two primitives.
     * @param[in] t The parameterization of the two primitives.
     */
    void compute_pDpT(VectorXd& pDpT, const VectorXd& s, const VectorXd& t) const;

    /**
     * @brief Compute the partial derivative of the distance w.r.t. s.
     * @param[out] pDpS \f$ \frac{\partial D}{\partial s} \f$
     * @param[in] s The state of the two primitives.
     * @param[in] t The parameterization of the two primitives.
     */
    void compute_pDpS(VectorXd& pDpS, const VectorXd& s, const VectorXd& t) const;

    /**
     * @brief Compute the second partial derivative of the distance w.r.t. t.
     * @param[out] p2DpT2 \f$ \frac{\partial^2 D}{\partial t^2} \f$
     * @param[in] s The state of the two primitives.
     * @param[in] t The parameterization of the two primitives.
     */
    void compute_p2DpT2(MatrixXd& p2DpT2, const VectorXd& s, const VectorXd& t) const;

    /**
     * @brief Compute the second partial derivative of the distance w.r.t. s.
     * @param[out] p2DpS2 \f$ \frac{\partial^2 D}{\partial s^2} \f$
     * @param[in] s The state of the two primitives.
     * @param[in] t The parameterization of the two primitives.
     */
    void compute_p2DpS2(MatrixXd& p2DpS2, const VectorXd& s, const VectorXd& t) const;

    /**
     * @brief Compute the mixed partial derivative of the distance w.r.t. t and s.
     * @param[out] p2DpTpS \f$ \frac{\partial^2 D}{\partial t \partial s} \f$
     * @param[in] s The state of the two primitives.
     * @param[in] t The parameterization of the two primitives.
     */
    void compute_p2DpTpS(MatrixXd& p2DpTpS, const VectorXd& s, const VectorXd& t) const;

    /** @} */  // end of group

    /**
     * @name O
     * @brief O and its derivatives
     */
    /** @{ */  // start of group

    /**
     * @brief Compute the objective value.
     * @param[in] s The state of the two primitives.
     * @param[in] t The parameterization of the two primitives.
     * @return The objective value.
     */
    double compute_O(const VectorXd& s, const VectorXd& t) const;

    /**
     * @brief Compute the partial derivative of the objective value w.r.t. t.
     * @param[out] pOpT \f$ \frac{\partial O}{\partial t} \f$
     * @param[in] s The state of the two primitives.
     * @param[in] t The parameterization of the two primitives.
     */
    void compute_pOpT(VectorXd& pOpT, const VectorXd& s, const VectorXd& t) const;

    /**
     * @brief Compute the partial derivative of the objective value w.r.t. s.
     * @param[out] pOpS \f$ \frac{\partial O}{\partial s} \f$
     * @param[in] s The state of the two primitives.
     * @param[in] t The parameterization of the two primitives.
     */
    void compute_pOpS(VectorXd& pOpS, const VectorXd& s, const VectorXd& t) const;

    /**
     * @brief Compute the second partial derivative of the objective value w.r.t. t.
     * @param[out] p2OpT2 \f$ \frac{\partial^2 O}{\partial t^2} \f$
     * @param[in] s The state of the two primitives.
     * @param[in] t The parameterization of the two primitives.
     */
    void compute_p2OpT2(MatrixXd& p2OpT2, const VectorXd& s, const VectorXd& t) const;

    /**
     * @brief Compute the second partial derivative of the objective value w.r.t. s.
     * @param[out] p2OpS2 \f$ \frac{\partial^2 O}{\partial s^2} \f$
     * @param[in] s The state of the two primitives.
     * @param[in] t The parameterization of the two primitives.
     */
    void compute_p2OpS2(MatrixXd& p2OpS2, const VectorXd& s, const VectorXd& t) const;

    /**
     * @brief Compute the mixed partial derivative of the objective value w.r.t. t and s.
     * @param[out] p2OpTpS \f$ \frac{\partial^2 O}{\partial t \partial s} \f$
     * @param[in] s The state of the two primitives.
     * @param[in] t The parameterization of the two primitives.
     */
    void compute_p2OpTpS(MatrixXd& p2OpTpS, const VectorXd& s, const VectorXd& t) const;

    /** @} */  // end of group

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    //--- Finite difference tests
    void test_pOpT_WithFD(const VectorXd& s, const VectorXd& t) const;
    void test_pOpS_WithFD(const VectorXd& s, const VectorXd& t) const;
    void test_p2OpT2_WithFD(const VectorXd& s, const VectorXd& t) const;
    void test_p2OpS2_WithFD(const VectorXd& s, const VectorXd& t) const;
    void test_p2OpTpS_WithFD(const VectorXd& s, const VectorXd& t) const;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    /**
     * @brief Return both sizes of t.
     * @return Both sizes of t.
     */
    std::pair<int, int> get_sizes_t() const;
    /**
     * @brief Return the sum of t sizes.
     * @return The sum of t sizes.
     */
    int get_total_size_t() const;

    /**
     * @brief Get a string of both primitive descriptions
     * @param[in] primitive_A The first primitive.
     * @param[in] primitive_B The second primitive.
     * @return The string.
     */
    static std::string get_primitives_description(primitive_t primitive_A, primitive_t primitive_B);

private:
    /**
     * @brief Helper for the first or second primitive.
     */
    enum PRIMITIVE { A, B };

    /**
     * @brief Get s for a primitive.
     * @param[in] s The stacked s vector.
     * @param[in] P the primitive to use.
     * @return The state belonging to the given primitive.
     */
    Vector6d get_s(const VectorXd& s, PRIMITIVE P) const;

    /**
     * @brief Get s for a primitive.
     * @param[in] t The stacked s vector.
     * @param[in] P the primitive to use.
     * @return The parameterization belonging to the given primitive.
     */
    VectorXd get_t(const VectorXd& t, PRIMITIVE P) const;

    /**
     * @brief Check s and t.
     * @param[in] s The stacked s vector.
     * @param[in] t the stacked t vector.
     * @throws std::runtime_error If the size is wrong.
     */
    void check_inputs(const VectorXd& s, const VectorXd& t) const;

public:
    double regularizerWeight = 0.01;  ///< Regularitaion wheight for the objective.
    double constraintWeight = 10.0;   ///< Constraint weight for the objective.

    SoftUnilateralConstraint unilateralConstraint = SoftUnilateralConstraint(0.0, 1.0, 1e-3);  ///< The soft unilateral constraint.

    const primitive_t primitive_A;  ///< The first primitive.
    const primitive_t primitive_B;  ///< The seonc primitive.
};

}  // namespace DCA
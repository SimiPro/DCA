#pragma once

#include <DCA/Opt/FiniteDifference.h>
#include <DCA/Utils/Utils.h>

#include <array>
#include <variant>
#include <vector>

namespace DCA {

/**
 * @brief Definition of a %Primitive
 */

class Primitive : public Opt::FiniteDifference {
public:
    Primitive(const std::string& description, const double& radius);
    virtual ~Primitive() {}

    //Compute point on primitive (and its derivatives) based on s and t
    virtual Vector3d compute_P(const Vector6d& s, const VectorXd& t) const = 0;

    virtual Eigen::Matrix<double, 3, 6> compute_dPdS(const Vector6d& s, const VectorXd& t) const = 0;
    virtual Eigen::Matrix<double, 3, -1> compute_dPdT(const Vector6d& s, const VectorXd& t) const = 0;

    virtual std::array<Eigen::Matrix<double, 3, 6>, 6> compute_d2PdS2(const Vector6d& s, const VectorXd& t) const = 0;
    virtual std::vector<Eigen::Matrix<double, 3, 6>> compute_d2PdSdT(const Vector6d& s, const VectorXd& t) const = 0;

    //Helpers
    virtual Vector6d get_s() const = 0;
    virtual int SIZE_T() const = 0;

    //Helpers for pair generator
    Vector3d get_center_point() const;
    virtual double get_largest_dimension_from_center() const = 0;

    //Test derivatives with FD
    void test_dPdS_WithFD() const;
    void test_dPdT_WithFD() const;
    void test_d2PdS2_WithFD() const;
    void test_d2PdSdT_WithFD() const;
    void test_derivatives() const;

protected:
    //Check size of t
    void check_t(const VectorXd& t) const;

private:
    //Finite difference
    typedef std::vector<MatrixXd> TensorXd;
    typedef std::function<void(TensorXd& tensor, const VectorXd& s, const VectorXd& t)> Function_Tensor;

    void estimateTensor(TensorXd& tensor, const Function_Matrix& evaluate, const VectorXd& s, const VectorXd& t, D_WRT d_wrt, int firstDim,
                        int secondDim) const;
    void testTensor(const Function_Matrix& evaluate, const Function_Tensor& analytic, const VectorXd& s, const VectorXd& t, D_WRT d_wrt,
                    const std::string& name, int firstDim, int secondDim) const;
    void printFDCheck(const TensorXd& estimate, const TensorXd& analytic, const std::string& name) const;

    VectorXd get_t_forFD() const;

public:
    double radius;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Definition of a %Sphere
 */
class Sphere : public Primitive {
public:
    Sphere(const Vector3d& position, const double& radius);
    ~Sphere() {}

    Vector3d compute_P(const Vector6d& s, const VectorXd& t) const override;

    Eigen::Matrix<double, 3, 6> compute_dPdS(const Vector6d& s, const VectorXd& t) const override;
    Eigen::Matrix<double, 3, -1> compute_dPdT(const Vector6d& s, const VectorXd& t) const override;

    std::array<Eigen::Matrix<double, 3, 6>, 6> compute_d2PdS2(const Vector6d& s, const VectorXd& t) const override;
    std::vector<Eigen::Matrix<double, 3, 6>> compute_d2PdSdT(const Vector6d& s, const VectorXd& t) const override;

    Vector6d get_s() const override;
    int SIZE_T() const override;
    double get_largest_dimension_from_center() const override;

public:
    Vector3d position;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Definition of a %Capsule
 */
class Capsule : public Primitive {
public:
    Capsule(const Vector3d& startPosition, const Vector3d& endPosition, const double& radius);
    ~Capsule() {}

    Vector3d compute_P(const Vector6d& s, const VectorXd& t) const override;

    Eigen::Matrix<double, 3, 6> compute_dPdS(const Vector6d& s, const VectorXd& t) const override;
    Eigen::Matrix<double, 3, -1> compute_dPdT(const Vector6d& s, const VectorXd& t) const override;

    std::array<Eigen::Matrix<double, 3, 6>, 6> compute_d2PdS2(const Vector6d& s, const VectorXd& t) const override;
    std::vector<Eigen::Matrix<double, 3, 6>> compute_d2PdSdT(const Vector6d& s, const VectorXd& t) const override;

    Vector6d get_s() const override;
    int SIZE_T() const override;
    double get_largest_dimension_from_center() const override;

public:
    Vector3d startPosition, endPosition;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Definition of a %Rectangle
 */
class Rectangle : public Primitive {
public:
    Rectangle(const Vector3d& center, const Matrix3d& orientation, const Vector2d& dimensions, const double& radius);
    ~Rectangle() {}

    Vector3d compute_P(const Vector6d& s, const VectorXd& t) const override;

    Eigen::Matrix<double, 3, 6> compute_dPdS(const Vector6d& s, const VectorXd& t) const override;
    Eigen::Matrix<double, 3, -1> compute_dPdT(const Vector6d& s, const VectorXd& t) const override;

    std::array<Eigen::Matrix<double, 3, 6>, 6> compute_d2PdS2(const Vector6d& s, const VectorXd& t) const override;
    std::vector<Eigen::Matrix<double, 3, 6>> compute_d2PdSdT(const Vector6d& s, const VectorXd& t) const override;

    Vector6d get_s() const override;
    int SIZE_T() const override;
    double get_largest_dimension_from_center() const override;

private:
    std::pair<Vector3d, Vector3d> getLocalVectors() const;

public:
    Vector3d center;
    Matrix3d orientation;
    Vector2d dimensions;  //Convention: local dimensions in x and z direction
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Definition of a %Box
 */
class Box : public Primitive {
public:
    Box(const Vector3d& center, const Matrix3d& orientation, const Vector3d& dimensions, const double& radius);
    ~Box() {}

    Vector3d compute_P(const Vector6d& s, const VectorXd& t) const override;

    Eigen::Matrix<double, 3, 6> compute_dPdS(const Vector6d& s, const VectorXd& t) const override;
    Eigen::Matrix<double, 3, -1> compute_dPdT(const Vector6d& s, const VectorXd& t) const override;

    std::array<Eigen::Matrix<double, 3, 6>, 6> compute_d2PdS2(const Vector6d& s, const VectorXd& t) const override;
    std::vector<Eigen::Matrix<double, 3, 6>> compute_d2PdSdT(const Vector6d& s, const VectorXd& t) const override;

    Vector6d get_s() const override;
    int SIZE_T() const override;
    double get_largest_dimension_from_center() const override;

private:
    std::tuple<Vector3d, Vector3d, Vector3d> getLocalVectors() const;

public:
    Vector3d center;
    Matrix3d orientation;
    Vector3d dimensions;  //Convention: local dimensions in x, y, z direction
};

/**
 * @brief All possible primitives
 */
using primitive_t = std::variant<Sphere, Capsule, Rectangle, Box>;

}  // namespace DCA
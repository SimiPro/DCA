#pragma once

#include <DCA/FiniteDifference.h>
#include <DCA/Utils.h>

#include <array>
#include <variant>
#include <vector>

namespace DCA {

/**
 * @brief Definition of a %Primitive
 * 
 * Each %Primitive can be described by a 6-dimensional state
 * and a n-dimensional parameterization vector.
 * 
 * This means, any convex primitive which can be described in a continuously differentiable manner can be implemented.
 * 
 * The inheritance on FiniteDifference is just to make the automatic derivative checking easier.
 */
class Primitive : public FiniteDifference {
public:
    /**
     * @brief Create a primitive.
     * 
     * @param[in] The description of the primitive.
     * @param[in] safetyMargin The safety margin of the primitive.
     * @throws std::runtime_error if safety margin < 0
     */
    Primitive(const std::string& description, const double& safetyMargin);

    /**
     * @brief Default deconstructor.
     */
    virtual ~Primitive() = default;

    /**
     * @brief Compute a point on the primitive.
     * 
     * Compute a point on the primitive based on s and t.
     * @param[in] s The state of the primitive.
     * @param[in] t The parameterization of the point.
     * @return The point.
     */
    virtual Vector3d compute_P(const Vector6d& s, const VectorXd& t) const = 0;

    /**
     * @brief Compute the derivative of a point on the primitive with respect to s.
     * 
     * Compute the derivative of a point on the primitive based on s and t.
     * @param[in] s The state of the primitive.
     * @param[in] t The parameterization of the point.
     * @return The derivative \f$ \frac{dP}{ds} \f$.
     */
    virtual Eigen::Matrix<double, 3, 6> compute_dPdS(const Vector6d& s, const VectorXd& t) const = 0;

    /**
     * @brief Compute the derivative of a point on the primitive with respect to t.
     * 
     * Compute the derivative of a point on the primitive based on s and t.
     * @param[in] s The state of the primitive.
     * @param[in] t The parameterization of the point.
     * @return The derivative \f$ \frac{dP}{dt} \f$.
     */
    virtual Eigen::Matrix<double, 3, -1> compute_dPdT(const Vector6d& s, const VectorXd& t) const = 0;

    /**
     * @brief Compute the second derivative of a point on the primitive with respect to s.
     * 
     * Compute the second derivative of a point on the primitive based on s and t.
     * @param[in] s The state of the primitive.
     * @param[in] t The parameterization of the point.
     * @return The second derivative \f$ \frac{d^2P}{ds^2} \f$ as an array (size 6).
     */
    virtual std::array<Eigen::Matrix<double, 3, 6>, 6> compute_d2PdS2(const Vector6d& s, const VectorXd& t) const = 0;

    /**
     * @brief Compute the mixed second derivative of a point on the primitive with respect to s and t.
     * 
     * Compute the mixed second derivative of a point on the primitive based on s and t.
     * @param[in] s The state of the primitive.
     * @param[in] t The parameterization of the point.
     * @return The derivative \f$ \frac{d^2P}{ds dt} \f$ as a vector (size t).
     */
    virtual std::vector<Eigen::Matrix<double, 3, 6>> compute_d2PdSdT(const Vector6d& s, const VectorXd& t) const = 0;

    /**
     * @brief Helper to get the state.
     * @return The state of the primitive.
     */
    virtual Vector6d get_s() const = 0;

    /**
     * @brief Helper to get the size of t.
     * @return The number of dimensions needed for parameterization.
     */
    virtual int SIZE_T() const = 0;

    /**
     * @brief Get the center point.
     * 
     * Helper for pair generation.
     * @return The center point of the primitive.
     */
    Vector3d get_center_point() const;

    /**
     * @brief Get the largest dimension.
     * 
     * Helper for pair generation.
     * @return The largest dimension from the center point.
     */
    virtual double get_largest_dimension_from_center() const = 0;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    // Test derivatives with FD.
    void test_dPdS_WithFD() const;
    void test_dPdT_WithFD() const;
    void test_d2PdS2_WithFD() const;
    void test_d2PdSdT_WithFD() const;
    void test_derivatives() const;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

protected:
    /**
     * @brief Check the size of t.
     * @param[in] t The t vector.
     * @throw std::runtime_error If the size is wrong.
     */
    void check_t(const VectorXd& t) const;

private:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    VectorXd get_t_forFD() const;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

public:
    double safetyMargin;  ///< Internal storage of the safety margin.
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Definition of a %Sphere
 * 
 * The state of a sphere is as follows:
 * (x, y, z, -, -, -), which means that it is the center of the sphere and the other three parameters are not used.
 * 
 * The parameterization of a sphere is as follows:
 * (-), which means, it does not need any (0-dimensional).
 */
class Sphere : public Primitive {
public:
    /**
     * Create a sphere primitive.
     * 
     * @param[in] position The position of the sphere.
     * @param[in] radius The safety margin around the center.
     */
    Sphere(const Vector3d& position, const double& radius);

    /**
     * @brief Default deconstructor.
     */
    ~Sphere() = default;

    /**
     * @copydoc Primitive::compute_P
     * 
     * The paramterization t is not used, since it has size 0!
     */
    Vector3d compute_P(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::compute_dPdS
     */
    Eigen::Matrix<double, 3, 6> compute_dPdS(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::compute_dPdT
     * 
     * In the case of a sphere, this has size 3x0.
     */
    Eigen::Matrix<double, 3, -1> compute_dPdT(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::compute_d2PdS2
     */
    std::array<Eigen::Matrix<double, 3, 6>, 6> compute_d2PdS2(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::compute_d2PdSdT
     * 
     * In the case of a sphere, this returns an empty vector.
     */
    std::vector<Eigen::Matrix<double, 3, 6>> compute_d2PdSdT(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::get_s
     */
    Vector6d get_s() const override;

    /**
     * @copydoc Primitive::SIZE_T
     * 
     * In the case of a sphere, this returns 0.
     */
    int SIZE_T() const override;

    /**
     * @copydoc Primitive::get_largest_dimension_from_center
     * 
     * In the case of a sphere, this returns the radius.
     */
    double get_largest_dimension_from_center() const override;

public:
    Vector3d position;  ///< The position of the sphere in 3d space.
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Definition of a %Capsule
 * 
 * The state of a capsule is as follows:
 * (x1, y1, z1, x2, y2, z2), which means that it is the start and end point of the underlying line segment, stacked.
 * 
 * The parameterization of a capsule is as follows:
 * (t1), which means, the position on the underlying line segment from -1 to 1.
 */
class Capsule : public Primitive {
public:
    /**
     * Create a capsule primitive.
     * 
     * @param[in] startPosition The start position of the underlying line segment.
     * @param[in] endPosition The end position of the underlying line segment.
     * @param[in] radius The safety margin around the center.
     * @throws std::runtime_error If the radius is smaller than the distance between start and end.
     */
    Capsule(const Vector3d& startPosition, const Vector3d& endPosition, const double& radius);

    /**
     * @brief Default deconstructor.
     */
    ~Capsule() = default;

    /**
     * @copydoc Primitive::compute_P
     */
    Vector3d compute_P(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::compute_dPdS
     */
    Eigen::Matrix<double, 3, 6> compute_dPdS(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::compute_dPdT
     * 
     * In the case of a capsule, this has size 3x1.
     */
    Eigen::Matrix<double, 3, -1> compute_dPdT(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::compute_d2PdS2
     */
    std::array<Eigen::Matrix<double, 3, 6>, 6> compute_d2PdS2(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::compute_d2PdSdT
     * 
     * In the case of a capsule, this returns a vector with one matrix.
     */
    std::vector<Eigen::Matrix<double, 3, 6>> compute_d2PdSdT(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::get_s
     */
    Vector6d get_s() const override;

    /**
     * @copydoc Primitive::SIZE_T
     * 
     * In the case of a capsule, this returns 1.
     */
    int SIZE_T() const override;

    /**
     * @copydoc Primitive::get_largest_dimension_from_center
     */
    double get_largest_dimension_from_center() const override;

public:
    Vector3d startPosition;  ///< The start position of the capsule in 3d space.
    Vector3d endPosition;    ///< The end position of the capsule in 3d space.
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Definition of a %Rectangle
 * 
 * The state of a rectangle is as follows:
 * (x1, y1, z1, th1, th2, th3), which means that it is center of mass and the exponential map theta.
 * 
 * The parameterization of a rectangle is as follows:
 * (t1, t2), which means, the parameterisation of a bounded plane.
 */
class Rectangle : public Primitive {
public:
    /**
     * Create a rectangle primitive.
     * 
     * @param[in] center The center of mass.
     * @param[in] orientation The orientation of the rectangle.
     * @param[in] dimensions The two dimensions of the rectangle in x and z direction.
     * @param[in] safetyMargin The safety margin of the rectangle.
     * @throws std::runtime_error If safety margin < 0.
     * @throws std::runtime_error If any dimension < 0.
     * @throws std::runtime_error If the orientation is invalid.
     */
    Rectangle(const Vector3d& center, const Matrix3d& orientation, const Vector2d& dimensions, const double& safetyMargin = 0.001);
    /**
     * @brief Default deconstructor.
     */
    ~Rectangle() = default;

    /**
     * @copydoc Primitive::compute_P
     */
    Vector3d compute_P(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::compute_dPdS
     */
    Eigen::Matrix<double, 3, 6> compute_dPdS(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::compute_dPdT
     * 
     * In the case of a rectangle, this has size 3x2.
     */
    Eigen::Matrix<double, 3, -1> compute_dPdT(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::compute_d2PdS2
     */
    std::array<Eigen::Matrix<double, 3, 6>, 6> compute_d2PdS2(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::compute_d2PdSdT
     * 
     * In the case of a rectangle, this returns a vector with size 2.
     */
    std::vector<Eigen::Matrix<double, 3, 6>> compute_d2PdSdT(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::get_s
     */
    Vector6d get_s() const override;

    /**
     * @copydoc Primitive::SIZE_T
     * 
     * In the case of a rectangle, this returns 2.
     */
    int SIZE_T() const override;

    /**
     * @copydoc Primitive::get_largest_dimension_from_center
     */
    double get_largest_dimension_from_center() const override;

private:
    /**
     * @brief Get vectors from the coordinates.
     * @return The local vectors.
     */
    std::pair<Vector3d, Vector3d> getLocalVectors() const;

public:
    Vector3d center;       ///< The center of mass of the rectangle.
    Matrix3d orientation;  ///< The orientation of the rectangle.
    Vector2d dimensions;   ///< Convention: local dimensions in x and z direction
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Definition of a %Box
 * 
 * The state of a box is as follows:
 * (x1, y1, z1, th1, th2, th3), which means that it is center of mass and the exponential map theta.
 * 
 * The parameterization of a box is as follows:
 * (t1, t2, t3), which means, the parameterisation of a box.
 */
class Box : public Primitive {
public:
    /**
     * Create a box primitive.
     * 
     * @param[in] center The center of mass.
     * @param[in] orientation The orientation of the box.
     * @param[in] dimensions The three dimensions of the box in x, y and z direction.
     * @param[in] safetyMargin The safety margin of the box.
     * @throws std::runtime_error If safety margin < 0.
     * @throws std::runtime_error If any dimension < 0.
     * @throws std::runtime_error If the orientation is invalid.
     */
    Box(const Vector3d& center, const Matrix3d& orientation, const Vector3d& dimensions, const double& safetyMargin = 0.001);

    /**
     * @brief Default deconstructor.
     */
    ~Box() = default;

    /**
     * @copydoc Primitive::compute_P
     */
    Vector3d compute_P(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::compute_dPdS
     */
    Eigen::Matrix<double, 3, 6> compute_dPdS(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::compute_dPdT
     * 
     * In the case of a box, this has size 3x3.
     */
    Eigen::Matrix<double, 3, -1> compute_dPdT(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::compute_d2PdS2
     */
    std::array<Eigen::Matrix<double, 3, 6>, 6> compute_d2PdS2(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::compute_d2PdSdT
     * 
     * In the case of a box, this returns a vector with size 3.
     */
    std::vector<Eigen::Matrix<double, 3, 6>> compute_d2PdSdT(const Vector6d& s, const VectorXd& t) const override;

    /**
     * @copydoc Primitive::get_s
     */
    Vector6d get_s() const override;

    /**
     * @copydoc Primitive::SIZE_T
     * 
     * In the case of a rectangle, this returns 3.
     */
    int SIZE_T() const override;

    /**
     * @copydoc Primitive::get_largest_dimension_from_center
     */
    double get_largest_dimension_from_center() const override;

private:
    /**
     * @brief Get vectors from the coordinates.
     * @return The local vectors.
     */
    std::tuple<Vector3d, Vector3d, Vector3d> getLocalVectors() const;

public:
    Vector3d center;       ///< The center of mass of the box.
    Matrix3d orientation;  ///< The orientation of the box.
    Vector3d dimensions;   ///< Convention: local dimensions in x, y, z direction
};

/**
 * @brief All possible primitives
 */
using primitive_t = std::variant<Sphere, Capsule, Rectangle, Box>;

}  // namespace DCA
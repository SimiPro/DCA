#include <DCA/Utils/ExpCoords.h>
#include <DCA/Utils/Primitives.h>

namespace DCA {

Vector3d Primitive::get_center_point() const {
    Vector6d s = get_s();
    VectorXd t = VectorXd::Zero(SIZE_T());
    return compute_P(s, t);
}

Primitive::Primitive(const std::string& description, const double& safetyMargin) : Opt::FiniteDifference(description), safetyMargin(safetyMargin) {
    if (safetyMargin < 0.0)
        throw std::runtime_error("Error in Primitive::Primitive -> invalid safety margin");
    throw std::runtime_error("Error in Primitive::Primitive -> invalid safety margin");
}

void Primitive::test_dPdS_WithFD() const {
    auto eval = [&](VectorXd& vec, const VectorXd& s, const VectorXd& t) -> void { vec = compute_P(s, t); };
    auto anal = [&](MatrixXd& mat, const VectorXd& s, const VectorXd& t) -> void { mat = compute_dPdS(s, t); };

    Vector6d s = get_s();
    VectorXd t = get_t_forFD();
    testMatrix(eval, anal, s, t, D_S, "dPdS", 3);
}

void Primitive::test_dPdT_WithFD() const {
    auto eval = [&](VectorXd& vec, const VectorXd& s, const VectorXd& t) -> void { vec = compute_P(s, t); };
    auto anal = [&](MatrixXd& mat, const VectorXd& s, const VectorXd& t) -> void { mat = compute_dPdT(s, t); };

    Vector6d s = get_s();
    VectorXd t = get_t_forFD();
    testMatrix(eval, anal, s, t, D_T, "dPdT", 3);
}

void Primitive::test_d2PdS2_WithFD() const {
    auto eval = [&](MatrixXd& mat, const VectorXd& s, const VectorXd& t) -> void { mat = compute_dPdS(s, t); };
    auto anal = [&](TensorXd& ten, const VectorXd& s, const VectorXd& t) -> void {
        ten.clear();
        std::array<Eigen::Matrix<double, 3, 6>, 6> d2PdS2 = compute_d2PdS2(s, t);
        for (const auto& mat : d2PdS2)
            ten.emplace_back(mat);
    };

    Vector6d s = get_s();
    VectorXd t = get_t_forFD();
    testTensor(eval, anal, s, t, D_S, "d2PdS2", 3, 6);
}

void Primitive::test_d2PdSdT_WithFD() const {
    auto eval = [&](MatrixXd& mat, const VectorXd& s, const VectorXd& t) -> void { mat = compute_dPdS(s, t); };
    auto anal = [&](TensorXd& ten, const VectorXd& s, const VectorXd& t) -> void {
        ten.clear();
        std::vector<Eigen::Matrix<double, 3, 6>> d2PdSdT = compute_d2PdSdT(s, t);
        for (const auto& mat : d2PdSdT)
            ten.emplace_back(mat);
    };

    Vector6d s = get_s();
    VectorXd t = get_t_forFD();
    testTensor(eval, anal, s, t, D_T, "d2PdSdT", 3, 6);
}

void Primitive::test_derivatives() const {
    test_dPdS_WithFD();
    test_dPdT_WithFD();
    test_d2PdSdT_WithFD();
    test_d2PdS2_WithFD();
}

void Primitive::check_t(const VectorXd& t) const {
    if (t.size() != SIZE_T())
        throw std::runtime_error("ERROR in Primitive: input t does not have the correct size!");
}

VectorXd Primitive::get_t_forFD() const {
    return VectorXd::Random(SIZE_T());
}

void Primitive::estimateTensor(TensorXd& tensor, const Function_Matrix& evaluate, const VectorXd& s, const VectorXd& t, D_WRT d_wrt, int firstDim,
                               int secondDim) const {
    VectorXd p_der = (d_wrt == D_S) ? s : t;
    tensor.clear();

    auto eval = [&](MatrixXd& mat, D_WRT d_wrt) -> void {
        if (d_wrt == D_S)
            evaluate(mat, p_der, t);
        else if (d_wrt == D_T)
            evaluate(mat, s, p_der);
    };

    MatrixXd f_P(firstDim, secondDim), f_M(firstDim, secondDim);
    for (int i = 0; i < p_der.size(); i++) {
        f_P.setZero();
        f_M.setZero();

        p_der[i] += deltaFD;
        eval(f_P, d_wrt);

        p_der[i] -= 2.0 * deltaFD;
        eval(f_M, d_wrt);

        p_der[i] += deltaFD;
        tensor.emplace_back((f_P - f_M) / (2 * deltaFD));
    }
}

void Primitive::testTensor(const Function_Matrix& evaluate, const Function_Tensor& analytic, const VectorXd& s, const VectorXd& t, D_WRT d_wrt,
                           const std::string& name, int firstDim, int secondDim) const {
    VectorXd p_der = (d_wrt == D_S) ? s : t;
    uint pSize = (uint)p_der.size();

    TensorXd fdTensor(pSize, MatrixXd::Zero(firstDim, secondDim));
    TensorXd anaTensor(pSize, MatrixXd::Zero(firstDim, secondDim));

    setFDCheckIsBeingApplied(true);
    estimateTensor(fdTensor, evaluate, s, t, d_wrt, firstDim, secondDim);
    analytic(anaTensor, s, t);
    setFDCheckIsBeingApplied(false);

    printFDCheck(fdTensor, anaTensor, name);
}

void Primitive::printFDCheck(const TensorXd& estimate, const TensorXd& analytic, const std::string& name) const {
    bool checkSuccessful = true;
    Logger::print(Logger::BLUE,
                  "------------------------------------------------------------"
                  "-----------------------------------------------\n");
    Logger::print(Logger::DEFAULT, "Test ");
    Logger::print(Logger::CYAN, "%s", name.c_str());
    Logger::print(Logger::DEFAULT, " of ");
    Logger::print(Logger::YELLOW, "%s", description.c_str());
    if (analytic.size() > 0 && estimate.size() > 0) {
        double norm_ana = 0.0, norm_est = 0.0;
        for (int k = 0; k < (int)estimate.size(); k++) {
            norm_ana += analytic[k].norm();
            norm_est += estimate[k].norm();
        }
        Logger::print(Logger::DEFAULT, " with FD -> norms: Ana: %lf, FD: %lf\n", norm_ana, norm_est);
    } else {
        Logger::print(Logger::DEFAULT, " with FD -> norms: Ana: %lf, FD: %lf\n", 0.0, 0.0);
    }

    for (int k = 0; k < (int)estimate.size(); k++)
        for (int i = 0; i < (int)estimate[k].rows(); i++) {
            for (int j = 0; j < (int)estimate[k].cols(); j++) {
                double absErr = std::abs(estimate[k].coeff(i, j) - analytic[k].coeff(i, j));
                double relError = 2 * absErr / (eps + std::abs(estimate[k].coeff(i, j)) + std::abs(analytic[k].coeff(i, j)));
                if (relError > relTol && absErr > absTol) {
                    checkSuccessful = false;
                    Logger::print(Logger::RED, "   MISMATCH");
                    Logger::print(Logger::DEFAULT,
                                  " -> Element: (%d, %d, %d), Ana val: %lf, FD val: "
                                  "%lf. Error: %lf(%lf%%)\n",
                                  i, j, k, analytic[k].coeff(i, j), estimate[k].coeff(i, j), absErr, relError * 100);
                }
            }
        }
    if (checkSuccessful)
        Logger::print(Logger::GREEN, "   PASSED\n");
    Logger::print(Logger::DEFAULT, "End of check of %s\n", description.c_str());
    Logger::print(Logger::BLUE,
                  "------------------------------------------------------------"
                  "-----------------------------------------------\n");
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------

Sphere::Sphere(const Vector3d& position, const double& radius) : Primitive("Sphere", radius), position(position) {}

Vector3d Sphere::compute_P(const Vector6d& s, const VectorXd& t) const {
    return s.segment<3>(0);
}

Eigen::Matrix<double, 3, 6> Sphere::compute_dPdS(const Vector6d& s, const VectorXd& t) const {
    Eigen::Matrix<double, 3, 6> dPdS;
    dPdS.block<3, 3>(0, 0).setIdentity();
    dPdS.block<3, 3>(0, 3).setZero();
    return dPdS;
}

Eigen::Matrix<double, 3, -1> Sphere::compute_dPdT(const Vector6d& s, const VectorXd& t) const {
    Eigen::Matrix<double, 3, 0> dPdT;
    dPdT.setZero();
    return dPdT;
}

Vector6d Sphere::get_s() const {
    Vector6d s = Vector6d::Zero();
    s.segment<3>(0) = position;
    return s;
}

std::vector<Eigen::Matrix<double, 3, 6>> Sphere::compute_d2PdSdT(const Vector6d& s, const VectorXd& t) const {
    return {};
}

std::array<Eigen::Matrix<double, 3, 6>, 6> Sphere::compute_d2PdS2(const Vector6d& s, const VectorXd& t) const {
    static const Eigen::Matrix<double, 3, 6> mat = Eigen::Matrix<double, 3, 6>::Zero();
    return {mat, mat, mat, mat, mat, mat};
}

int Sphere::SIZE_T() const {
    return 0;
}

double Sphere::get_largest_dimension_from_center() const {
    return safetyMargin;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------

Capsule::Capsule(const Vector3d& startPosition, const Vector3d& endPosition, const double& radius)
    : Primitive("Capsule", radius), startPosition(startPosition), endPosition(endPosition) {
    if (radius > (startPosition - endPosition).norm())
        throw std::runtime_error("Error in Capsule::Capsule -> invalid radius");
}

Vector3d Capsule::compute_P(const Vector6d& s, const VectorXd& t) const {
    check_t(t);
    return s.segment<3>(0) + 0.5 * (t[0] + 1.0) * (s.segment<3>(3) - s.segment<3>(0));
}

Eigen::Matrix<double, 3, 6> Capsule::compute_dPdS(const Vector6d& s, const VectorXd& t) const {
    check_t(t);
    Eigen::Matrix<double, 3, 6> dPdS;
    dPdS.block<3, 3>(0, 0) = (1.0 - 0.5 * (t[0] + 1.0)) * Matrix3d::Identity();
    dPdS.block<3, 3>(0, 3) = 0.5 * (t[0] + 1.0) * Matrix3d::Identity();
    return dPdS;
}

Eigen::Matrix<double, 3, -1> Capsule::compute_dPdT(const Vector6d& s, const VectorXd& t) const {
    Eigen::Matrix<double, 3, 1> dPdT = 0.5 * (s.segment<3>(3) - s.segment<3>(0));
    return dPdT;
}

std::array<Eigen::Matrix<double, 3, 6>, 6> Capsule::compute_d2PdS2(const Vector6d& s, const VectorXd& t) const {
    static const Eigen::Matrix<double, 3, 6> mat = Eigen::Matrix<double, 3, 6>::Zero();
    return {mat, mat, mat, mat, mat, mat};
}

std::vector<Eigen::Matrix<double, 3, 6>> Capsule::compute_d2PdSdT(const Vector6d& s, const VectorXd& t) const {
    Eigen::Matrix<double, 3, 6> mat;
    mat.block<3, 3>(0, 0) = -0.5 * Matrix3d::Identity();
    mat.block<3, 3>(0, 3) = 0.5 * Matrix3d::Identity();
    return {mat};
}

Vector6d Capsule::get_s() const {
    Vector6d s;
    s << startPosition, endPosition;
    return s;
}

int Capsule::SIZE_T() const {
    return 1;
}

double Capsule::get_largest_dimension_from_center() const {
    return 0.5 * (startPosition - endPosition).norm() + safetyMargin;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------

Rectangle::Rectangle(const Vector3d& center, const Matrix3d& orientation, const Vector2d& dimensions, const double& safetyMargin)
    : Primitive("Rectangle", safetyMargin), center(center), orientation(orientation), dimensions(dimensions) {
    if (fabs(orientation.determinant() - 1.0) > 1e-6)
        throw std::runtime_error("Error in Rectangle::Rectangle -> invalid orientation");

    if (dimensions[0] < 0.0 || dimensions[1] < 0.0)
        throw std::runtime_error("Error in Rectangle::Rectangle -> invalid dimensions");
}

std::pair<Vector3d, Vector3d> Rectangle::getLocalVectors() const {
    return {Vector3d::UnitX() * dimensions[0] / 2.0, Vector3d::UnitZ() * dimensions[1] / 2.0};
}

Vector3d Rectangle::compute_P(const Vector6d& s, const VectorXd& t) const {
    check_t(t);
    Vector3d com = s.segment<3>(0);
    Vector3d theta = s.segment<3>(3);  // this exp. coord. describe orientation of the box in the world

    // local dims of our cube
    auto [u1L, u2L] = getLocalVectors();

    // now transform them into world coords
    Vector3d u1W = ExpCoords::get_w(theta, u1L);
    Vector3d u2W = ExpCoords::get_w(theta, u2L);

    return com + t[0] * u1W + t[1] * u2W;
}

Eigen::Matrix<double, 3, 6> Rectangle::compute_dPdS(const Vector6d& s, const VectorXd& t) const {
    Eigen::Matrix<double, 3, 6> dpds;

    dpds.setZero();
    Vector3d com = s.segment<3>(0);
    Vector3d theta = s.segment<3>(3);  // this exp. coord. describe orientation of the box in the world

    // local dims of our cube
    auto [u1L, u2L] = getLocalVectors();

    // now transform them into world coords
    Matrix3d du1W_ds = ExpCoords::get_dwdr(theta, u1L);
    Matrix3d du2W_ds = ExpCoords::get_dwdr(theta, u2L);

    dpds.block<3, 3>(0, 0).setIdentity();
    dpds.block<3, 3>(0, 3) = t[0] * du1W_ds + t[1] * du2W_ds;
    return dpds;
}

Eigen::Matrix<double, 3, -1> Rectangle::compute_dPdT(const Vector6d& s, const VectorXd& t) const {
    check_t(t);
    Vector3d com = s.segment<3>(0);
    Vector3d theta = s.segment<3>(3);  // this exp. coord. describe orientation of the box in the world

    // local dims of our cube
    auto [u1L, u2L] = getLocalVectors();

    // now transform them into world coords
    Vector3d u1W = ExpCoords::get_w(theta, u1L);
    Vector3d u2W = ExpCoords::get_w(theta, u2L);

    Eigen::Matrix<double, 3, 2> res;
    res.col(0) = u1W;
    res.col(1) = u2W;

    return res;
}

std::array<Eigen::Matrix<double, 3, 6>, 6> Rectangle::compute_d2PdS2(const Vector6d& s, const VectorXd& t) const {
    static const Eigen::Matrix<double, 3, 6> mat = Eigen::Matrix<double, 3, 6>::Zero();
    std::array<Eigen::Matrix<double, 3, 6>, 6> d2PdS2 = {mat, mat, mat, mat, mat, mat};

    Vector3d theta = s.segment<3>(3);  // this exp. coord. describe orientation of the box in the world

    // local dims of our cube
    auto [u1L, u2L] = getLocalVectors();

    // now transform them into world coords
    {
        Matrix3d ddu1W_ds2 = ExpCoords::get_ddwdr_dr1(theta, u1L);
        Matrix3d ddu2W_ds2 = ExpCoords::get_ddwdr_dr1(theta, u2L);
        d2PdS2[3].block<3, 3>(0, 3) = t[0] * ddu1W_ds2 + t[1] * ddu2W_ds2;
    }
    {
        Matrix3d ddu1W_ds2 = ExpCoords::get_ddwdr_dr2(theta, u1L);
        Matrix3d ddu2W_ds2 = ExpCoords::get_ddwdr_dr2(theta, u2L);
        d2PdS2[4].block<3, 3>(0, 3) = t[0] * ddu1W_ds2 + t[1] * ddu2W_ds2;
    }
    {
        Matrix3d ddu1W_ds2 = ExpCoords::get_ddwdr_dr3(theta, u1L);
        Matrix3d ddu2W_ds2 = ExpCoords::get_ddwdr_dr3(theta, u2L);
        d2PdS2[5].block<3, 3>(0, 3) = t[0] * ddu1W_ds2 + t[1] * ddu2W_ds2;
    }

    return d2PdS2;
}

std::vector<Eigen::Matrix<double, 3, 6>> Rectangle::compute_d2PdSdT(const Vector6d& s, const VectorXd& t) const {
    Eigen::Matrix<double, 3, 6> dpds;

    dpds.setZero();
    Vector3d com = s.segment<3>(0);
    Vector3d theta = s.segment<3>(3);  // this exp. coord. describe orientation of the rectangle in the world

    // local dims of our cube
    auto [u1L, u2L] = getLocalVectors();

    // now transform them into world coords
    Matrix3d du1W_ds = ExpCoords::get_dwdr(theta, u1L);
    Matrix3d du2W_ds = ExpCoords::get_dwdr(theta, u2L);

    dpds.block<3, 3>(0, 0).setZero();
    Eigen::Matrix<double, 3, 6> dpds1 = dpds, dpds2 = dpds;
    dpds1.block<3, 3>(0, 3) = du1W_ds;
    dpds2.block<3, 3>(0, 3) = du2W_ds;
    return {dpds1, dpds2};
}

Vector6d Rectangle::get_s() const {
    Vector6d s = Vector6d::Zero();
    s.segment<3>(0) = center;
    s.segment<3>(3) = ExpCoords::get_theta(orientation);
    return s;
}

int Rectangle::SIZE_T() const {
    return 2;
}

double Rectangle::get_largest_dimension_from_center() const {
    return 0.5 * dimensions.norm();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------

Box::Box(const Vector3d& center, const Matrix3d& orientation, const Vector3d& dimensions, const double& safetyMargin)
    : Primitive("Box", safetyMargin), center(center), orientation(orientation), dimensions(dimensions) {
    if (fabs(orientation.determinant() - 1.0) > 1e-6)
        throw std::runtime_error("Error in Box::Box -> invalid orientation");

    if (dimensions[0] < 0.0 || dimensions[1] < 0.0 || dimensions[2] < 0.0)
        throw std::runtime_error("Error in Box::Box -> invalid dimensions");
}

std::tuple<Vector3d, Vector3d, Vector3d> Box::getLocalVectors() const {
    return {Vector3d::UnitX() * dimensions[0] / 2.0, Vector3d::UnitY() * dimensions[1] / 2.0, Vector3d::UnitZ() * dimensions[2] / 2.0};
}

Vector3d Box::compute_P(const Vector6d& s, const VectorXd& t) const {
    check_t(t);
    Vector3d com = s.segment<3>(0);
    Vector3d theta = s.segment<3>(3);  // this exp. coord. describe orientation of the box in the world

    // local dims of our cube
    auto [u1L, u2L, u3L] = getLocalVectors();

    // now transform them into world coords
    Vector3d u1W = ExpCoords::get_w(theta, u1L);
    Vector3d u2W = ExpCoords::get_w(theta, u2L);
    Vector3d u3W = ExpCoords::get_w(theta, u3L);

    return com + t[0] * u1W + t[1] * u2W + t[2] * u3W;
}

Eigen::Matrix<double, 3, 6> Box::compute_dPdS(const Vector6d& s, const VectorXd& t) const {
    Eigen::Matrix<double, 3, 6> dpds;

    dpds.setZero();
    Vector3d com = s.segment<3>(0);
    Vector3d theta = s.segment<3>(3);  // this exp. coord. describe orientation of the box in the world

    // local dims of our cube
    auto [u1L, u2L, u3L] = getLocalVectors();

    // now transform them into world coords
    Matrix3d du1W_ds = ExpCoords::get_dwdr(theta, u1L);
    Matrix3d du2W_ds = ExpCoords::get_dwdr(theta, u2L);
    Matrix3d du3W_ds = ExpCoords::get_dwdr(theta, u3L);

    dpds.block<3, 3>(0, 0).setIdentity();
    dpds.block<3, 3>(0, 3) = t[0] * du1W_ds + t[1] * du2W_ds + t[2] * du3W_ds;
    return dpds;
}

Eigen::Matrix<double, 3, -1> Box::compute_dPdT(const Vector6d& s, const VectorXd& t) const {
    check_t(t);
    Vector3d com = s.segment<3>(0);
    Vector3d theta = s.segment<3>(3);  // this exp. coord. describe orientation of the box in the world

    // local dims of our cube
    auto [u1L, u2L, u3L] = getLocalVectors();

    // now transform them into world coords
    Vector3d u1W = ExpCoords::get_w(theta, u1L);
    Vector3d u2W = ExpCoords::get_w(theta, u2L);
    Vector3d u3W = ExpCoords::get_w(theta, u3L);

    Eigen::Matrix<double, 3, 3> res;
    res.col(0) = u1W;
    res.col(1) = u2W;
    res.col(2) = u3W;

    return res;
}

std::array<Eigen::Matrix<double, 3, 6>, 6> Box::compute_d2PdS2(const Vector6d& s, const VectorXd& t) const {
    static const Eigen::Matrix<double, 3, 6> mat = Eigen::Matrix<double, 3, 6>::Zero();
    std::array<Eigen::Matrix<double, 3, 6>, 6> d2PdS2 = {mat, mat, mat, mat, mat, mat};

    Vector3d theta = s.segment<3>(3);  // this exp. coord. describe orientation of the box in the world

    // local dims of our cube
    auto [u1L, u2L, u3L] = getLocalVectors();

    // now transform them into world coords
    {
        Matrix3d ddu1W_ds2 = ExpCoords::get_ddwdr_dr1(theta, u1L);
        Matrix3d ddu2W_ds2 = ExpCoords::get_ddwdr_dr1(theta, u2L);
        Matrix3d ddu3W_ds2 = ExpCoords::get_ddwdr_dr1(theta, u3L);
        d2PdS2[3].block<3, 3>(0, 3) = t[0] * ddu1W_ds2 + t[1] * ddu2W_ds2 + t[2] * ddu3W_ds2;
    }
    {
        Matrix3d ddu1W_ds2 = ExpCoords::get_ddwdr_dr2(theta, u1L);
        Matrix3d ddu2W_ds2 = ExpCoords::get_ddwdr_dr2(theta, u2L);
        Matrix3d ddu3W_ds2 = ExpCoords::get_ddwdr_dr2(theta, u3L);
        d2PdS2[4].block<3, 3>(0, 3) = t[0] * ddu1W_ds2 + t[1] * ddu2W_ds2 + t[2] * ddu3W_ds2;
    }
    {
        Matrix3d ddu1W_ds2 = ExpCoords::get_ddwdr_dr3(theta, u1L);
        Matrix3d ddu2W_ds2 = ExpCoords::get_ddwdr_dr3(theta, u2L);
        Matrix3d ddu3W_ds2 = ExpCoords::get_ddwdr_dr3(theta, u3L);
        d2PdS2[5].block<3, 3>(0, 3) = t[0] * ddu1W_ds2 + t[1] * ddu2W_ds2 + t[2] * ddu3W_ds2;
    }

    return d2PdS2;
}

std::vector<Eigen::Matrix<double, 3, 6>> Box::compute_d2PdSdT(const Vector6d& s, const VectorXd& t) const {
    Eigen::Matrix<double, 3, 6> dpds;

    dpds.setZero();
    Vector3d com = s.segment<3>(0);
    Vector3d theta = s.segment<3>(3);  // this exp. coord. describe orientation of the box in the world

    // local dims of our cube
    auto [u1L, u2L, u3L] = getLocalVectors();

    // now transform them into world coords
    Matrix3d du1W_ds = ExpCoords::get_dwdr(theta, u1L);
    Matrix3d du2W_ds = ExpCoords::get_dwdr(theta, u2L);
    Matrix3d du3W_ds = ExpCoords::get_dwdr(theta, u3L);

    dpds.block<3, 3>(0, 0).setZero();
    Eigen::Matrix<double, 3, 6> dpds1 = dpds, dpds2 = dpds, dpds3 = dpds;
    dpds1.block<3, 3>(0, 3) = du1W_ds;
    dpds2.block<3, 3>(0, 3) = du2W_ds;
    dpds3.block<3, 3>(0, 3) = du3W_ds;
    return {dpds1, dpds2, dpds3};
}

Vector6d Box::get_s() const {
    Vector6d s = Vector6d::Zero();
    s.segment<3>(0) = center;
    s.segment<3>(3) = ExpCoords::get_theta(orientation);
    return s;
}

int Box::SIZE_T() const {
    return 3;
}

double Box::get_largest_dimension_from_center() const {
    return 0.5 * dimensions.norm();
}
}  // namespace DCA
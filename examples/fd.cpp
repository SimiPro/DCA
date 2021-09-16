#include <DCA/Interactions/Solver.h>

int main(int argc, char const *argv[]) {
    using namespace DCA;

    auto frand = [&]() -> double { return (double)rand() / RAND_MAX; };
    auto rrand = [&]() -> Matrix3d {
        return Eigen::AngleAxisd(frand(), Vector3d::UnitX()).matrix() * Eigen::AngleAxisd(frand(), Vector3d::UnitY()).matrix() *
               Eigen::AngleAxisd(frand(), Vector3d::UnitZ()).matrix();
    };
    auto dim2rand = [&]() -> Vector2d { return Vector2d(fabs(frand()), fabs(frand())); };
    auto dim3rand = [&]() -> Vector3d { return Vector3d(fabs(frand()), fabs(frand()), fabs(frand())); };

    Sphere(Vector3d::Random(), 0.0).test_derivatives();
    Capsule(Vector3d::Random(), Vector3d::Random(), 0.0).test_derivatives();
    Rectangle(Vector3d::Random(), rrand(), dim2rand(), 0.0).test_derivatives();
    Box(Vector3d::Random(), rrand(), dim3rand(), 0.0).test_derivatives();

    Interactions::Solver(Sphere(Vector3d::Random(), 0.0), Sphere(Vector3d::Random(), 0.0)).test_derivatives();
    Interactions::Solver(Sphere(Vector3d::Random(), 0.0), Capsule(Vector3d::Random(), Vector3d::Random(), 0.0)).test_derivatives();
    Interactions::Solver(Sphere(Vector3d::Random(), 0.0), Rectangle(Vector3d::Random(), rrand(), dim2rand(), 0.0)).test_derivatives();
    Interactions::Solver(Sphere(Vector3d::Random(), 0.0), Box(Vector3d::Random(), rrand(), dim3rand(), 0.0)).test_derivatives();

    Interactions::Solver(Capsule(Vector3d::Random(), Vector3d::Random(), 0.0), Capsule(Vector3d::Random(), Vector3d::Random(), 0.0)).test_derivatives();
    Interactions::Solver(Capsule(Vector3d::Random(), Vector3d::Random(), 0.0), Rectangle(Vector3d::Random(), rrand(), dim2rand(), 0.0)).test_derivatives();
    Interactions::Solver(Capsule(Vector3d::Random(), Vector3d::Random(), 0.0), Box(Vector3d::Random(), rrand(), dim3rand(), 0.0)).test_derivatives();

    Interactions::Solver(Rectangle(Vector3d::Random(), rrand(), dim2rand(), 0.0), Rectangle(Vector3d::Random(), rrand(), dim2rand(), 0.0)).test_derivatives();
    Interactions::Solver(Rectangle(Vector3d::Random(), rrand(), dim2rand(), 0.0), Box(Vector3d::Random(), rrand(), dim3rand(), 0.0)).test_derivatives();

    Interactions::Solver(Box(Vector3d::Random(), rrand(), dim3rand(), 0.0), Box(Vector3d::Random(), rrand(), dim3rand(), 0.0)).test_derivatives();
}

// Tell the library that we want to run FD checks.
// This will automatically make functions public,
// and set solver residuals.
#define RUN_FD_CHECK

#include <DCA/CapsuleVsCapsule.h>
#include <DCA/SphereVsCapsule.h>
#include <DCA/SphereVsSphere.h>

int main(int argc, char const *argv[]) {
    using namespace DCA;

    Vector6d P6;
    P6 << 1, 4, 45, 56, 1, 3;
    Vector9d P9;
    P9 << P6, 4, 78, 2;
    Vector12d P12;
    P12 << P9, 4, 7, 6;

    Vector2d props;  // radius for both primitives
    props << 1, 0.8;

    CapsuleDistanceObjective obj;
    Vector2d X;
    CapsuleVsCapsule::solveForX(P12, X);

    // Checking pairwise
    SphereVsSphere::check_dDdP_3(P6, props);
    SphereVsSphere::check_dDdP_6(P6, props);

    SphereVsSphere::check_d2DdP2_3(P6, props);
    SphereVsSphere::check_d2DdP2_6(P6, props);

    SphereVsCapsule::check_dDdP_3(P9, props);
    SphereVsCapsule::check_dDdP_6(P9, props);
    SphereVsCapsule::check_dDdP_9(P9, props);

    SphereVsCapsule::check_d2DdP2_3(P9, props);
    SphereVsCapsule::check_d2DdP2_6(P9, props);
    SphereVsCapsule::check_d2DdP2_9(P9, props);

    CapsuleVsCapsule::check_dDdP_6(P12, props);
    CapsuleVsCapsule::check_dDdP_12(P12, props);

    // The Hessians are only approximations, therefore the FD check fails.
    CapsuleVsCapsule::check_d2DdP2_6(P12, props);
    CapsuleVsCapsule::check_d2DdP2_12(P12, props);

    CapsuleVsCapsule::check_dXdP_12(P12, X);

    // Checking CapsuleDistanceObjective
    obj.check_dOdX_2(X, P12);
    obj.check_d2OdX2_2(X, P12);

    obj.check_dDdX_2(X, P12);
    obj.check_d2DdX2_2(X, P12);

    obj.check_dDdP_12(P12, X);
    obj.check_d2DdP2_12(P12, X);

    obj.check_d2DdXdP_12(P12, X);

    return 0;
}

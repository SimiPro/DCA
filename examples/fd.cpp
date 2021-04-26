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

    Vector2d props;
    props << 1, 0.8;

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

    CapsuleVsCapsule::check_d2DdP2_6(P12, props);
    CapsuleVsCapsule::check_d2DdP2_12(P12, props);

    return 0;
}

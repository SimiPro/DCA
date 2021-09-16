#pragma once

#include <DCA/Pair.h>
#include <DCA/Utils/Primitives.h>

namespace DCA {

/**
 * @brief Public %API.
 * 
 * This is the main public %API which should be used.
 */
class API {
public:
    static double compute_D(const primitive_t& p_a, const primitive_t& p_b);
    static void compute_dDdS(Vector12d& dDdS, const primitive_t& p_a, const primitive_t& p_b);
    static void compute_d2DdS2(Matrix12d& d2DdS2, const primitive_t& p_a, const primitive_t& p_b);

    static double compute_D(const primitive_t& p_a, const primitive_t& p_b, const VectorXd& t);
    static void compute_dDdS(Vector12d& dDdS, const primitive_t& p_a, const primitive_t& p_b, const VectorXd& t);
    static void compute_d2DdS2(Matrix12d& d2DdS2, const primitive_t& p_a, const primitive_t& p_b, const VectorXd& t);

    static std::pair<Vector3d, Vector3d> compute_closest_points(const primitive_t& p_a, const primitive_t& p_b);
    static std::pair<Vector3d, Vector3d> compute_closest_points(const primitive_t& p_a, const primitive_t& p_b, const VectorXd& t);

    static void compute_t(VectorXd& t, const primitive_t& p_a, const primitive_t& p_b);
};

}  // namespace DCA
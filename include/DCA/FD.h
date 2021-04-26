#ifndef __DCA__FD_H__
#define __DCA__FD_H__

#include <iostream>

#include "utils.h"

#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_DEFAULT "\x1b[0m"

namespace DCA {

#define FD_CHECK_dDdP(NumDer, NumDofs, NumProps, StartDer, NAME)             \
    static void check_dDdP_##NumDer(                                         \
        const Eigen::Matrix<double, NumDofs, 1>& P,                          \
        const Eigen::Matrix<double, NumProps, 1>& props) {                   \
        FD_Check<NumDer, NumDofs, NumProps, StartDer>::check_dDdP(           \
            P, props,                                                        \
            [&](const Eigen::Matrix<double, NumDofs, 1>& P,                  \
                const Eigen::Matrix<double, NumProps, 1>& props) -> double { \
                return compute_D(P, props);                                  \
            },                                                               \
            [&](Eigen::Matrix<double, NumDer, 1>& dDdP,                      \
                const Eigen::Matrix<double, NumDofs, 1>& P,                  \
                const Eigen::Matrix<double, NumProps, 1>& props) -> void {   \
                compute_dDdP(dDdP, P, props);                                \
            },                                                               \
            NAME);                                                           \
    }

#define FD_CHECK_d2DdP2(NumDer, NumDofs, NumProps, StartDer, NAME)         \
    static void check_d2DdP2_##NumDer(                                     \
        const Eigen::Matrix<double, NumDofs, 1>& P,                        \
        const Eigen::Matrix<double, NumProps, 1>& props) {                 \
        FD_Check<NumDer, NumDofs, NumProps, StartDer>::check_d2DdP2(       \
            P, props,                                                      \
            [&](Eigen::Matrix<double, NumDer, 1>& dDdP,                    \
                const Eigen::Matrix<double, NumDofs, 1>& P,                \
                const Eigen::Matrix<double, NumProps, 1>& props) -> void { \
                compute_dDdP(dDdP, P, props);                              \
            },                                                             \
            [&](Eigen::Matrix<double, NumDer, NumDer>& d2DdP2,             \
                const Eigen::Matrix<double, NumDofs, 1>& P,                \
                const Eigen::Matrix<double, NumProps, 1>& props) -> void { \
                compute_d2DdP2(d2DdP2, P, props);                          \
            },                                                             \
            NAME);                                                         \
    }

template <int NumDer, int NumDofs, int NumProps, int StartDer>
class FD_Check {
public:
    using dof_v = Eigen::Matrix<double, NumDofs, 1>;
    using props_v = Eigen::Matrix<double, NumProps, 1>;
    using der_v = Eigen::Matrix<double, NumDer, 1>;
    using der_m = Eigen::Matrix<double, NumDer, NumDer>;

    using evaluate_f =
        std::function<double(const dof_v& P, const props_v& props)>;
    using dDdP_f =
        std::function<void(der_v& dDdP, const dof_v& P, const props_v& props)>;

    using evaluate_fv =
        std::function<void(der_v& der, const dof_v& P, const props_v& props)>;
    using d2DdP2_fv = std::function<void(der_m& d2DdP2, const dof_v& P,
                                         const props_v& props)>;

    static void estimate_dDdP(der_v& dDdP, const dof_v& P, const props_v& props,
                              evaluate_f evaluate) {
        dof_v p_der(P);

        double f_P, f_M;
        for (uint i = StartDer; i < StartDer + NumDer; i++) {
            double tmpVal = p_der(i);
            p_der(i) = tmpVal + deltaFD;

            f_P = evaluate(p_der, props);
            p_der(i) = tmpVal - deltaFD;

            f_M = evaluate(p_der, props);
            p_der(i) = tmpVal;

            dDdP(i - StartDer) = (f_P - f_M) / (2 * deltaFD);
        }
    }

    static void check_dDdP(const dof_v& P, const props_v& props,
                           evaluate_f evaluate, dDdP_f analytic,
                           const char* name) {
        der_v analytic_dDdP;
        analytic(analytic_dDdP, P, props);

        der_v fd_dDdP;
        estimate_dDdP(fd_dDdP, P, props, evaluate);

        bool hasMismatch = false;
        std::cout << ANSI_COLOR_CYAN << "Checking " << name
                  << ANSI_COLOR_DEFAULT << std::endl;
        for (int i = 0; i < NumDer; i++) {
            if (fabs(analytic_dDdP(i) - fd_dDdP(i)) > 1e-5) {
                std::cout << ANSI_COLOR_RED << i
                          << " does not match. (Analytic: " << analytic_dDdP(i)
                          << ", FD: " << fd_dDdP(i) << ")" << ANSI_COLOR_DEFAULT
                          << std::endl;
                hasMismatch = true;
            }
        }

        if (!hasMismatch) {
            std::cout << ANSI_COLOR_GREEN << "All Good." << ANSI_COLOR_DEFAULT
                      << std::endl;
        }
    }

    static void estimate_d2DdP2(der_m& d2DdP2, const dof_v& P,
                                const props_v& props, evaluate_fv evaluate) {
        dof_v p_der(P);

        der_v f_P, f_M;
        for (uint i = StartDer; i < NumDer + StartDer; i++) {
            double tmpVal = p_der(i);
            p_der(i) = tmpVal + deltaFD;

            evaluate(f_P, p_der, props);
            p_der(i) = tmpVal - deltaFD;

            evaluate(f_M, p_der, props);
            p_der(i) = tmpVal;

            d2DdP2.col(i - StartDer) = (f_P - f_M) / (2 * deltaFD);
        }
    }

    static void check_d2DdP2(const dof_v& P, const props_v& props,
                             evaluate_fv evaluate, d2DdP2_fv analytic,
                             const char* name) {
        der_m analytic_d2DdP2;
        analytic(analytic_d2DdP2, P, props);

        der_m fd_d2DdP2;
        estimate_d2DdP2(fd_d2DdP2, P, props, evaluate);

        bool hasMismatch = false;
        std::cout << ANSI_COLOR_CYAN << "Checking " << name
                  << ANSI_COLOR_DEFAULT << std::endl;
        for (int i = 0; i < NumDer; i++) {
            for (int j = 0; j < NumDer; j++) {
                if (fabs(analytic_d2DdP2(i, j) - fd_d2DdP2(i, j)) > 1e-5) {
                    std::cout << ANSI_COLOR_RED << i << "/" << j
                              << " does not match. (Analytic: "
                              << analytic_d2DdP2(i, j)
                              << ", FD: " << fd_d2DdP2(i, j) << ")"
                              << ANSI_COLOR_DEFAULT << std::endl;
                    hasMismatch = true;
                }
            }
        }

        if (!hasMismatch) {
            std::cout << ANSI_COLOR_GREEN << "All Good." << ANSI_COLOR_DEFAULT
                      << std::endl;
        }
    }
};

}  // namespace DCA

#endif
#ifndef __DCA__FD_H__
#define __DCA__FD_H__

#include <iostream>

#include "utils.h"

#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_DEFAULT "\x1b[0m"

namespace DCA {

/**
 * This helper is for the Newton Objective,
 * which is non static and uses different
 * function names (compute_O instead of compute_D)
 */
#define FD_CHECK_dOdX(NumDer, NumDofs, NumProps, StartDer, NAME)            \
    void check_dOdX_##NumDer(const Eigen::Matrix<double, NumDofs, 1>& X,    \
                             const Eigen::Matrix<double, NumProps, 1>& P) { \
        FD_Check<NumDer, NumDofs, NumProps, StartDer>::check_dDdP(          \
            X, P,                                                           \
            [&](const Eigen::Matrix<double, NumDofs, 1>& X,                 \
                const Eigen::Matrix<double, NumProps, 1>& P) -> double {    \
                return compute_O(P, X);                                     \
            },                                                              \
            [&](Eigen::Matrix<double, NumDer, 1>& dOdP,                     \
                const Eigen::Matrix<double, NumDofs, 1>& X,                 \
                const Eigen::Matrix<double, NumProps, 1>& P) -> void {      \
                compute_dOdX(dOdP, P, X);                                   \
            },                                                              \
            NAME);                                                          \
    }

#define FD_CHECK_d2OdX2(NumDer, NumDofs, NumProps, StartDer, NAME)            \
    void check_d2OdX2_##NumDer(const Eigen::Matrix<double, NumDofs, 1>& X,    \
                               const Eigen::Matrix<double, NumProps, 1>& P) { \
        FD_Check<NumDer, NumDofs, NumProps, StartDer>::check_d2DdP2(          \
            X, P,                                                             \
            [&](Eigen::Matrix<double, NumDer, 1>& dOdX,                       \
                const Eigen::Matrix<double, NumDofs, 1>& X,                   \
                const Eigen::Matrix<double, NumProps, 1>& P) -> void {        \
                return compute_dOdX(dOdX, P, X);                              \
            },                                                                \
            [&](Eigen::Matrix<double, NumDer, NumDer>& d2OdP2,                \
                const Eigen::Matrix<double, NumDofs, 1>& X,                   \
                const Eigen::Matrix<double, NumProps, 1>& P) -> void {        \
                compute_d2OdX2(d2OdP2, P, X);                                 \
            },                                                                \
            NAME);                                                            \
    }

#define FD_CHECK_dDdX(NumDer, NumDofs, NumProps, StartDer, NAME)            \
    void check_dDdX_##NumDer(const Eigen::Matrix<double, NumDofs, 1>& X,    \
                             const Eigen::Matrix<double, NumProps, 1>& P) { \
        FD_Check<NumDer, NumDofs, NumProps, StartDer>::check_dDdP(          \
            X, P,                                                           \
            [&](const Eigen::Matrix<double, NumDofs, 1>& X,                 \
                const Eigen::Matrix<double, NumProps, 1>& P) -> double {    \
                return compute_D(P, X);                                     \
            },                                                              \
            [&](Eigen::Matrix<double, NumDer, 1>& dOdP,                     \
                const Eigen::Matrix<double, NumDofs, 1>& X,                 \
                const Eigen::Matrix<double, NumProps, 1>& P) -> void {      \
                compute_dDdX(dOdP, P, X);                                   \
            },                                                              \
            NAME);                                                          \
    }

#define FD_CHECK_d2DdX2(NumDer, NumDofs, NumProps, StartDer, NAME)            \
    void check_d2DdX2_##NumDer(const Eigen::Matrix<double, NumDofs, 1>& X,    \
                               const Eigen::Matrix<double, NumProps, 1>& P) { \
        FD_Check<NumDer, NumDofs, NumProps, StartDer>::check_d2DdP2(          \
            X, P,                                                             \
            [&](Eigen::Matrix<double, NumDer, 1>& dOdX,                       \
                const Eigen::Matrix<double, NumDofs, 1>& X,                   \
                const Eigen::Matrix<double, NumProps, 1>& P) -> void {        \
                return compute_dDdX(dOdX, P, X);                              \
            },                                                                \
            [&](Eigen::Matrix<double, NumDer, NumDer>& d2OdP2,                \
                const Eigen::Matrix<double, NumDofs, 1>& X,                   \
                const Eigen::Matrix<double, NumProps, 1>& P) -> void {        \
                compute_d2DdX2(d2OdP2, P, X);                                 \
            },                                                                \
            NAME);                                                            \
    }

#define FD_CHECK_d2DdXdP(NumDer, NumDofs, NumProps, StartDer, NAME)            \
    void check_d2DdXdP_##NumDer(const Eigen::Matrix<double, NumDofs, 1>& P,    \
                                const Eigen::Matrix<double, NumProps, 1>& X) { \
        FD_Check<NumDer, NumDofs, NumProps, StartDer>::check_d2DdXdP(          \
            P, X,                                                              \
            [&](Eigen::Matrix<double, NumProps, 1>& dDdX,                      \
                const Eigen::Matrix<double, NumDofs, 1>& P,                    \
                const Eigen::Matrix<double, NumProps, 1>& X) -> void {         \
                compute_dDdX(dDdX, P, X);                                      \
            },                                                                 \
            [&](Eigen::Matrix<double, NumProps, NumDer>& d2DdXdP,              \
                const Eigen::Matrix<double, NumDofs, 1>& P,                    \
                const Eigen::Matrix<double, NumProps, 1>& X) -> void {         \
                compute_d2DdXdP(d2DdXdP, P, X);                                \
            },                                                                 \
            NAME);                                                             \
    }

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

#define FD_CHECK_dDdP_NON_STATIC(NumDer, NumDofs, NumProps, StartDer, NAME)  \
    void check_dDdP_##NumDer(                                                \
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

#define FD_CHECK_d2DdP2_NON_STATIC(NumDer, NumDofs, NumProps, StartDer, NAME) \
    void check_d2DdP2_##NumDer(                                               \
        const Eigen::Matrix<double, NumDofs, 1>& P,                           \
        const Eigen::Matrix<double, NumProps, 1>& props) {                    \
        FD_Check<NumDer, NumDofs, NumProps, StartDer>::check_d2DdP2(          \
            P, props,                                                         \
            [&](Eigen::Matrix<double, NumDer, 1>& dDdP,                       \
                const Eigen::Matrix<double, NumDofs, 1>& P,                   \
                const Eigen::Matrix<double, NumProps, 1>& props) -> void {    \
                compute_dDdP(dDdP, P, props);                                 \
            },                                                                \
            [&](Eigen::Matrix<double, NumDer, NumDer>& d2DdP2,                \
                const Eigen::Matrix<double, NumDofs, 1>& P,                   \
                const Eigen::Matrix<double, NumProps, 1>& props) -> void {    \
                compute_d2DdP2(d2DdP2, P, props);                             \
            },                                                                \
            NAME);                                                            \
    }

#define FD_CHECK_dXdP(NumDer, NumDofs, NumProps, StartDer, NAME)       \
    static void check_dXdP_##NumDer(                                   \
        const Eigen::Matrix<double, NumDofs, 1>& P,                    \
        const Eigen::Matrix<double, NumProps, 1>& X) {                 \
        FD_Check<NumDer, NumDofs, NumProps, StartDer>::check_dXdP(     \
            P, X,                                                      \
            [&](const Eigen::Matrix<double, NumDofs, 1>& P,            \
                Eigen::Matrix<double, NumProps, 1>& X) -> void {       \
                solveForX(P, X);                                       \
            },                                                         \
            [&](Eigen::Matrix<double, NumProps, NumDer>& dXdP,         \
                const Eigen::Matrix<double, NumDofs, 1>& P,            \
                const Eigen::Matrix<double, NumProps, 1>& X) -> void { \
                compute_dXdP(dXdP, P, X);                              \
            },                                                         \
            NAME);                                                     \
    }
template <int NumDer, int NumDofs, int NumProps, int StartDer>
class FD_Check {
public:
    using dof_v = Eigen::Matrix<double, NumDofs, 1>;
    using props_v = Eigen::Matrix<double, NumProps, 1>;
    using der_v = Eigen::Matrix<double, NumDer, 1>;
    using der_m = Eigen::Matrix<double, NumDer, NumDer>;
    using d2DdXdP_m = Eigen::Matrix<double, NumProps, NumDofs>;

    using evaluate_f =
        std::function<double(const dof_v& P, const props_v& props)>;
    using dDdP_f =
        std::function<void(der_v& dDdP, const dof_v& P, const props_v& props)>;

    using evaluate_fv =
        std::function<void(der_v& der, const dof_v& P, const props_v& props)>;
    using d2DdP2_fv = std::function<void(der_m& d2DdP2, const dof_v& P,
                                         const props_v& props)>;

    using evaluate_dXdP_fv = std::function<void(const dof_v& P, props_v& X)>;
    using evaluate_d2DdXdP_fv =
        std::function<void(props_v& der, const dof_v& P, const props_v& X)>;
    using d2DdXdP_fv =
        std::function<void(d2DdXdP_m&, const dof_v& P, const props_v& X)>;

    static double deltaFD() {
        static double deltaFD_ = 1e-6;
        return deltaFD_;
    }

    static void estimate_dDdP(der_v& dDdP, const dof_v& P, const props_v& props,
                              evaluate_f evaluate) {
        dof_v p_der(P);

        double f_P, f_M;
        for (int i = StartDer; i < StartDer + NumDer; i++) {
            double tmpVal = p_der(i);
            p_der(i) = tmpVal + deltaFD();

            f_P = evaluate(p_der, props);
            p_der(i) = tmpVal - deltaFD();

            f_M = evaluate(p_der, props);
            p_der(i) = tmpVal;

            dDdP(i - StartDer) = (f_P - f_M) / (2 * deltaFD());
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
        for (int i = StartDer; i < NumDer + StartDer; i++) {
            double tmpVal = p_der(i);
            p_der(i) = tmpVal + deltaFD();

            evaluate(f_P, p_der, props);
            p_der(i) = tmpVal - deltaFD();

            evaluate(f_M, p_der, props);
            p_der(i) = tmpVal;

            d2DdP2.col(i - StartDer) = (f_P - f_M) / (2 * deltaFD());
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

    static void estimate_d2DdXdP(d2DdXdP_m& d2DdXdP, const dof_v& P,
                                 const props_v& X,
                                 evaluate_d2DdXdP_fv evaluate) {
        dof_v p_der(P);

        props_v f_P, f_M;
        for (int i = StartDer; i < NumDer + StartDer; i++) {
            double tmpVal = p_der(i);
            p_der(i) = tmpVal + deltaFD();

            evaluate(f_P, p_der, X);
            p_der(i) = tmpVal - deltaFD();

            evaluate(f_M, p_der, X);
            p_der(i) = tmpVal;

            d2DdXdP.col(i - StartDer) = (f_P - f_M) / (2 * deltaFD());
        }
    }

    static void check_d2DdXdP(const dof_v& P, const props_v& X,
                              evaluate_d2DdXdP_fv evaluate, d2DdXdP_fv analytic,
                              const char* name) {
        d2DdXdP_m analytic_d2DdXdP;
        analytic(analytic_d2DdXdP, P, X);

        d2DdXdP_m fd_d2DdXdP;
        estimate_d2DdXdP(fd_d2DdXdP, P, X, evaluate);

        bool hasMismatch = false;
        std::cout << ANSI_COLOR_CYAN << "Checking " << name
                  << ANSI_COLOR_DEFAULT << std::endl;
        for (int i = 0; i < NumProps; i++) {
            for (int j = 0; j < NumDer; j++) {
                if (fabs(analytic_d2DdXdP(i, j) - fd_d2DdXdP(i, j)) > 1e-5) {
                    std::cout << ANSI_COLOR_RED << i << "/" << j
                              << " does not match. (Analytic: "
                              << analytic_d2DdXdP(i, j)
                              << ", FD: " << fd_d2DdXdP(i, j) << ")"
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

    static void estimate_dXdP(d2DdXdP_m& dXdP, const dof_v& P, const props_v& X,
                              evaluate_dXdP_fv evaluate) {
        dof_v p_der(P);

        props_v f_P, f_M;
        for (int i = StartDer; i < StartDer + NumDer; i++) {
            double tmpVal = p_der(i);
            p_der(i) = tmpVal + deltaFD();

            evaluate(p_der, f_P);
            p_der(i) = tmpVal - deltaFD();

            evaluate(p_der, f_M);
            p_der(i) = tmpVal;

            dXdP.col(i - StartDer) = (f_P - f_M) / (2 * deltaFD());
        }
    }

    static void check_dXdP(const dof_v& P, const props_v& X,
                           evaluate_dXdP_fv evaluate, d2DdXdP_fv analytic,
                           const char* name) {
        d2DdXdP_m analytic_dXdP;
        analytic(analytic_dXdP, P, X);

        d2DdXdP_m fd_dXdP;
        estimate_dXdP(fd_dXdP, P, X, evaluate);

        bool hasMismatch = false;
        std::cout << ANSI_COLOR_CYAN << "Checking " << name
                  << ANSI_COLOR_DEFAULT << std::endl;
        for (int i = 0; i < NumProps; i++) {
            for (int j = 0; j < NumDer; j++) {
                if (fabs(analytic_dXdP(i, j) - fd_dXdP(i, j)) > 1e-5) {
                    std::cout
                        << ANSI_COLOR_RED << i << "/" << j
                        << " does not match. (Analytic: " << analytic_dXdP(i, j)
                        << ", FD: " << fd_dXdP(i, j) << ")"
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
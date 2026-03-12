#include "spdlog/spdlog.h"
#include <cstdlib>

#include "gropt_utils.hpp"

#include "gropt_params.hpp"
#include "op_safe.hpp"

namespace Gropt {

static Eigen::VectorXd get_SAFE_compute(const Eigen::VectorXd &G, int Naxis, double dt, Op_SAFE &opF) {
    int N = static_cast<int>(G.size()) / Naxis;

    opF.init();

    Eigen::VectorXd temp;
    temp.setZero(opF.Ax_size);
    Eigen::VectorXd G_copy = G;
    opF.forward(G_copy, temp);
    opF.check(temp);

    opF.x_temp.setZero();
    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            if (opF.n_terms == 3) {
                opF.x_temp(j * N + i) = temp(j * opF.n_terms * N + i) + temp(j * opF.n_terms * N + i + N) +
                                        temp(j * opF.n_terms * N + i + 2 * N);
            } else if (opF.n_terms == 2) {
                opF.x_temp(j * N + i) = temp(j * opF.n_terms * N + i) + temp(j * opF.n_terms * N + i + N);
            }
        }
    }

    return opF.x_temp;
}

static ProblemData make_safe_pdata(const Eigen::VectorXd &G, int Naxis, double dt) {
    int N = static_cast<int>(G.size()) / Naxis;
    ProblemData pdata;
    pdata.dt = dt;
    pdata.N = N;
    pdata.Naxis = Naxis;
    pdata.inv_vec.setOnes(N * Naxis);
    pdata.set_vals.setZero(N * Naxis);
    pdata.set_vals.array() *= NAN;
    pdata.set_vals(0) = 0.0;
    pdata.set_vals(N - 1) = 0.0;
    pdata.fixer.setOnes(N * Naxis);
    pdata.fixer(0) = 0.0;
    pdata.fixer(N - 1) = 0.0;
    return pdata;
}

Eigen::VectorXd get_SAFE_eigen(const Eigen::VectorXd &G, int Naxis, double dt, bool true_safe, int new_first_axis) {
    ProblemData pdata = make_safe_pdata(G, Naxis, dt);
    Op_SAFE opF(pdata, 1.0, 1.0);
    opF.safe_params.set_demo_params();
    opF.safe_params.swap_first_axes(new_first_axis);
    return get_SAFE_compute(G, Naxis, dt, opF);
}

Eigen::VectorXd get_SAFE_eigen(const Eigen::VectorXd &G, int Naxis, double dt, bool true_safe, int new_first_axis,
                               const Eigen::VectorXd &tau1, const Eigen::VectorXd &tau2, const Eigen::VectorXd &tau3,
                               const Eigen::VectorXd &a1, const Eigen::VectorXd &a2, const Eigen::VectorXd &a3,
                               const Eigen::VectorXd &stim_limit, const Eigen::VectorXd &g_scale) {
    ProblemData pdata = make_safe_pdata(G, Naxis, dt);
    Op_SAFE opF(pdata, 1.0, 1.0);
    opF.safe_params.set_params(tau1, tau2, tau3, a1, a2, a3, stim_limit, g_scale);
    opF.safe_params.swap_first_axes(new_first_axis);
    return get_SAFE_compute(G, Naxis, dt, opF);
}

void test_eigen_assertions(int test_type) {
    switch (test_type) {
    case 1: {
        spdlog::info("Starting test_eigen_assertions(1)...");
        Eigen::VectorXd v(3);
        double x = v(10); // out of bounds
        spdlog::info("Finsihed test_eigen_assertions(1)  x = {}", x);
        break;
    }
    case 2: {
        spdlog::info("Starting test_eigen_assertions(2)...");
        Eigen::VectorXd v(5);
        if (!v.array().isNaN().any()) {
            spdlog::info("No NaN found");
            throw std::runtime_error("expected NaN initialization");
        } else {
            spdlog::info("NaN found");
        }
        spdlog::info("Finished test_eigen_assertions(2)  v[0] = {}", v[0]);
        break;
    }
    case 3: {
        spdlog::info("Starting test_eigen_assertions(3)...");
        Eigen::VectorXd a(3), b(5);
        Eigen::VectorXd c = a + b; // size mismatch
        spdlog::info("Finished test_eigen_assertions(3)  c[4] = {}", c[4]);
        break;
    }
    default:
        throw std::invalid_argument("unknown test type");
    }
}

} // namespace Gropt

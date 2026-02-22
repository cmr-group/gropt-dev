#include <cstdlib>
#include "spdlog/spdlog.h"

#include "gropt_utils.hpp"

#include "gropt_params.hpp"
#include "op_safe.hpp"

namespace Gropt {

void set_verbose(int level) {
    if (level == 0) {
        spdlog::set_level(spdlog::level::off);
    } else if (level == 1) {
        spdlog::set_level(spdlog::level::critical);
    } else if (level == 2) {
        spdlog::set_level(spdlog::level::err);
    } else if (level == 3) {
        spdlog::set_level(spdlog::level::warn);
    } else if (level == 4) {
        spdlog::set_level(spdlog::level::info);
    } else if (level == 5) {
        spdlog::set_level(spdlog::level::debug);
    } else {
        spdlog::set_level(spdlog::level::trace);
    }
}

void get_SAFE(int N, int Naxis, double dt, double *G_in,
              bool true_safe, int new_first_axis, bool demo_params,
              double *tau1, double *tau2, double *tau3,
              double *a1, double *a2, double *a3,
              double *stim_limit, double *g_scale,
              double **out, int &out_size)
{
    spdlog::trace("get_SAFE(): start");

    Eigen::VectorXd G;
    G.setZero(N);
    for (int i=0; i<N; i++) {
        G(i) = G_in[i];
    }

    spdlog::trace("get_SAFE(): copied G");

    ProblemData pdata;
    pdata.dt = dt;
    pdata.N = N;
    pdata.Naxis = Naxis;
    pdata.inv_vec.setOnes(N * Naxis);
    pdata.set_vals.setZero(N * Naxis);
    pdata.set_vals.array() *= NAN;
    pdata.set_vals(0) = 0.0;
    pdata.set_vals(N-1) = 0.0;
    pdata.fixer.setOnes(N * Naxis);
    pdata.fixer(0) = 0.0;
    pdata.fixer(N-1) = 0.0;

    spdlog::trace("get_SAFE(): finished pdata");

    Op_SAFE opF(pdata, 1.0, 1.0);
    if (demo_params) {
        opF.safe_params.set_demo_params();
    } else {
        opF.safe_params.set_params(tau1, tau2, tau3, a1, a2, a3, stim_limit, g_scale);
    }
    opF.safe_params.swap_first_axes(new_first_axis);
    opF.init();

    spdlog::trace("get_SAFE(): finished Op_SAFE");

    Eigen::VectorXd temp;
    temp.setZero(opF.Ax_size);
    opF.forward(G, temp);

    opF.check(temp);

    spdlog::trace("get_SAFE(): finished forward");

    opF.x_temp.setZero();
    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            opF.x_temp(j*N+i) = temp(j*3*N+i) + temp(j*3*N+i+N) + temp(j*3*N+i+2*N);
        }
    }

    spdlog::trace("get_SAFE(): finished combine");

    out_size = opF.x_temp.size();
    *out = (double*)std::malloc(out_size * sizeof(double));
    for (int i = 0; i < out_size; i++) {
        (*out)[i] = opF.x_temp(i);
    }

    spdlog::trace("get_SAFE(): finished copying out");

}

static Eigen::VectorXd get_SAFE_impl(const Eigen::VectorXd &G, int Naxis, double dt,
                                     bool true_safe, int new_first_axis, bool demo_params,
                                     double *tau1, double *tau2, double *tau3,
                                     double *a1, double *a2, double *a3,
                                     double *stim_limit, double *g_scale)
{
    int N = static_cast<int>(G.size()) / Naxis;

    ProblemData pdata;
    pdata.dt = dt;
    pdata.N = N;
    pdata.Naxis = Naxis;
    pdata.inv_vec.setOnes(N * Naxis);
    pdata.set_vals.setZero(N * Naxis);
    pdata.set_vals.array() *= NAN;
    pdata.set_vals(0) = 0.0;
    pdata.set_vals(N-1) = 0.0;
    pdata.fixer.setOnes(N * Naxis);
    pdata.fixer(0) = 0.0;
    pdata.fixer(N-1) = 0.0;

    Op_SAFE opF(pdata, 1.0, 1.0);
    if (demo_params) {
        opF.safe_params.set_demo_params();
    } else {
        opF.safe_params.set_params(tau1, tau2, tau3, a1, a2, a3, stim_limit, g_scale);
    }
    opF.safe_params.swap_first_axes(new_first_axis);
    opF.init();

    Eigen::VectorXd temp;
    temp.setZero(opF.Ax_size);
    Eigen::VectorXd G_copy = G;
    opF.forward(G_copy, temp);
    opF.check(temp);

    opF.x_temp.setZero();
    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            opF.x_temp(j*N+i) = temp(j*3*N+i) + temp(j*3*N+i+N) + temp(j*3*N+i+2*N);
        }
    }

    return opF.x_temp;
}

Eigen::VectorXd get_SAFE_eigen(const Eigen::VectorXd &G, int Naxis, double dt,
                               bool true_safe, int new_first_axis)
{
    return get_SAFE_impl(G, Naxis, dt, true_safe, new_first_axis, true,
                         nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
}

Eigen::VectorXd get_SAFE_eigen(const Eigen::VectorXd &G, int Naxis, double dt,
                               bool true_safe, int new_first_axis,
                               const Eigen::VectorXd &tau1, const Eigen::VectorXd &tau2, const Eigen::VectorXd &tau3,
                               const Eigen::VectorXd &a1, const Eigen::VectorXd &a2, const Eigen::VectorXd &a3,
                               const Eigen::VectorXd &stim_limit, const Eigen::VectorXd &g_scale)
{
    return get_SAFE_impl(G, Naxis, dt, true_safe, new_first_axis, false,
                         const_cast<double*>(tau1.data()), const_cast<double*>(tau2.data()), const_cast<double*>(tau3.data()),
                         const_cast<double*>(a1.data()), const_cast<double*>(a2.data()), const_cast<double*>(a3.data()),
                         const_cast<double*>(stim_limit.data()), const_cast<double*>(g_scale.data()));
}

} // namespace Gropt

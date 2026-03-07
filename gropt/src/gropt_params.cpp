#include "spdlog/spdlog.h"

#include "gropt_params.hpp"

#include "op_bvalue.hpp"
#include "op_gradient.hpp"
#include "op_identity.hpp"
#include "op_moment.hpp"
#include "op_safe.hpp"
#include "op_slew.hpp"
#include "op_tv.hpp"

namespace Gropt {

GroptParams::GroptParams() {}

void GroptParams::vec_init_simple(int _N, int _Naxis, double first_val, double last_val) {
    if (_N > 0) {
        N = _N;
    }

    if (_Naxis > 0) {
        Naxis = _Naxis;
    }

    Ntot = N * Naxis;

    pdata.inv_vec.setOnes(N * Naxis);

    pdata.set_vals.setZero(N * Naxis);
    pdata.set_vals.array() *= NAN;
    pdata.set_vals(0) = first_val;
    pdata.set_vals(N - 1) = last_val;

    pdata.fixer.setOnes(N * Naxis);
    pdata.fixer(0) = 0.0;
    pdata.fixer(N - 1) = 0.0;

    pdata.X0.setOnes(N * Naxis);
    pdata.X0 *= .01;
    pdata.X0(0) = first_val;
    pdata.X0(N - 1) = last_val;

    vec_init_status = N;
}

void GroptParams::setvec_X0(int _N, int _Naxis, double *_X0, bool set_others) {
    if (_N > 0) {
        N = _N;
    }

    if (_Naxis > 0) {
        Naxis = _Naxis;
    }

    Ntot = N * Naxis;

    pdata.X0.setOnes(N * Naxis);
    for (int i = 0; i < N * Naxis; i++) {
        pdata.X0(i) = _X0[i];
    }

    if (set_others) {
        pdata.inv_vec.setOnes(N * Naxis);

        pdata.set_vals.setZero(N * Naxis);
        pdata.set_vals.array() *= NAN;

        pdata.fixer.setOnes(N * Naxis);

        for (int j = 0; j < Naxis; j++) {
            pdata.set_vals((j * N)) = pdata.X0((j * N));
            pdata.set_vals((j * N) + N - 1) = pdata.X0((j * N) + N - 1);

            pdata.fixer((j * N)) = 0.0;
            pdata.fixer((j * N) + N - 1) = 0.0;
        }
    }

    vec_init_status = N;
}

void GroptParams::setvec_X0(const Eigen::VectorXd &_X0, int _Naxis, bool set_others) {
    int _N = static_cast<int>(_X0.size()) / _Naxis;
    // const_cast is safe here — the raw-pointer overload only reads from the data
    setvec_X0(_N, _Naxis, const_cast<double *>(_X0.data()), set_others);
}

void GroptParams::diff_init(double _dt, double _TE, double _T_90, double _T_180, double _T_readout) {
    dt = _dt;
    Naxis = 1;

    double T_90 = _T_90;
    double T_180 = _T_180;
    double T_readout = _T_readout;
    double TE = _TE;

    N = (int)((TE - T_readout) / dt) + 1;
    Ntot = N * Naxis;

    int ind_inv = (int)(TE / 2.0 / dt);
    pdata.inv_vec.setOnes(N);
    for (int i = ind_inv; i < N; i++) {
        pdata.inv_vec(i) = -1.0;
    }

    int ind_90_end, ind_180_start, ind_180_end;
    ind_90_end = ceil(T_90 / dt);
    ind_180_start = floor((TE / 2.0 - T_180 / 2.0) / dt);
    ind_180_end = ceil((TE / 2.0 + T_180 / 2.0) / dt);

    pdata.set_vals.setOnes(N);
    pdata.set_vals.array() *= NAN;
    for (int i = 0; i <= ind_90_end; i++) {
        pdata.set_vals(i) = 0.0;
    }
    for (int i = ind_180_start; i <= ind_180_end; i++) {
        pdata.set_vals(i) = 0.0;
    }
    pdata.set_vals(0) = 0.0;
    pdata.set_vals(N - 1) = 0.0;

    pdata.fixer.setOnes(N);
    for (int i = 0; i < pdata.set_vals.size(); i++) {
        if (!isnan(pdata.set_vals(i))) {
            pdata.fixer(i) = 0.0;
        }
    }

    pdata.X0.setOnes(N);
    for (int i = 0; i < pdata.set_vals.size(); i++) {
        if (!isnan(pdata.set_vals(i))) {
            pdata.X0(i) = pdata.set_vals(i);
        } else {
            pdata.X0(i) = 1e-2; // Initial value for non-fixed points
        }
    }

    vec_init_status = N;
}

void GroptParams::set_ils_solver(std::string _ils_method) {
    spdlog::info("set_ils_solver: {}", _ils_method);
    if (_ils_method == "CG") {
        ils_method = CG;
    } else if (_ils_method == "NLCG") {
        ils_method = NLCG;
    } else if (_ils_method == "BiCGstabl") {
        ils_method = BiCGstabl;
    } else {
        spdlog::error("Unknown Indirect Linear Solver method: {}", _ils_method);
        throw std::invalid_argument("Unknown Indirect Linear Solver method");
    }
}

void GroptParams::vec_reduce_simple(int N_reduce) {
    N -= N_reduce;
    Ntot = N * Naxis;

    pdata.inv_vec.setOnes(N * Naxis);

    Eigen::VectorXd set_vals_new;
    set_vals_new.setZero(N * Naxis);
    set_vals_new.array() *= NAN;
    set_vals_new(0) = pdata.set_vals(0);
    set_vals_new(N - 1) = pdata.set_vals(pdata.set_vals.size() - 1);

    Eigen::VectorXd fixer_new;
    fixer_new.setOnes(N * Naxis);
    fixer_new(0) = 0.0;
    fixer_new(N - 1) = 0.0;

    pdata.set_vals = set_vals_new;
    pdata.fixer = fixer_new;

    vec_init_status = N;
}

void GroptParams::prepare() {

    spdlog::trace("GroptParams::prepare() start");

    if (N * Naxis != Ntot) {
        Ntot = N * Naxis;
        spdlog::warn("Ntot is not consistent with N and Naxis");
        spdlog::warn("Setting Ntot = {}", Ntot);
    }

    if (vec_init_status != N) {
        spdlog::info("set_vals and inv_vec were not initialized, calling vec_init_simple()");
        vec_init_simple(-1, -1, 0.0, 0.0);
    }

    for (int i = 0; i < all_op.size(); i++) {
        all_op[i]->init();
    }
    for (int i = 0; i < all_obj.size(); i++) {
        all_obj[i]->init();
    }

    op_prep_status = N;

    spdlog::trace("GroptParams::prepare() end");
}

void GroptParams::add_gmax(double gmax, bool rot_variant, double weight_mod) {
    all_op.push_back(std::make_unique<Op_Gradient>(pdata, gmax, rot_variant, weight_mod));
}

void GroptParams::add_smax(double smax, bool rot_variant, double weight_mod) {
    all_op.push_back(std::make_unique<Op_Slew>(pdata, smax, rot_variant, weight_mod));
}

void GroptParams::add_moment(double order, double target, double tol0, std::string units, int moment_axis,
                             int start_idx0, int stop_idx0, int ref_idx0, double weight_mod) {
    all_op.push_back(std::make_unique<Op_Moment>(pdata, order, target, tol0, units, moment_axis, start_idx0, stop_idx0,
                                                 ref_idx0, weight_mod));
}

void GroptParams::add_SAFE(double stim_thresh, double *tau1, double *tau2, double *tau3, double *a1, double *a2,
                           double *a3, double *stim_limit, double *g_scale, int new_first_axis, bool demo_params,
                           double weight_mod) {
    auto op_F = std::make_unique<Op_SAFE>(pdata, stim_thresh, weight_mod);
    if (demo_params) {
        op_F->safe_params.set_demo_params();
    } else {
        op_F->safe_params.set_params(tau1, tau2, tau3, a1, a2, a3, stim_limit, g_scale);
    }
    op_F->safe_params.swap_first_axes(new_first_axis);
    all_op.push_back(std::move(op_F));
}

void GroptParams::add_SAFE_vec(int N_vec, double *stim_thresh_vec, double *tau1, double *tau2, double *tau3, double *a1,
                               double *a2, double *a3, double *stim_limit, double *g_scale, int new_first_axis,
                               bool demo_params, double weight_mod) {
    auto op_F = std::make_unique<Op_SAFE>(pdata, N_vec, stim_thresh_vec, weight_mod);
    if (demo_params) {
        op_F->safe_params.set_demo_params();
    } else {
        op_F->safe_params.set_params(tau1, tau2, tau3, a1, a2, a3, stim_limit, g_scale);
    }
    op_F->safe_params.swap_first_axes(new_first_axis);
    all_op.push_back(std::move(op_F));
}

void GroptParams::add_SAFE(double stim_thresh, int new_first_axis, double weight_mod) {
    auto op_F = std::make_unique<Op_SAFE>(pdata, stim_thresh, weight_mod);
    op_F->safe_params.set_demo_params();
    op_F->safe_params.swap_first_axes(new_first_axis);
    all_op.push_back(std::move(op_F));
}

void GroptParams::add_SAFE(double stim_thresh, const Eigen::VectorXd &tau1, const Eigen::VectorXd &tau2,
                           const Eigen::VectorXd &tau3, const Eigen::VectorXd &a1, const Eigen::VectorXd &a2,
                           const Eigen::VectorXd &a3, const Eigen::VectorXd &stim_limit, const Eigen::VectorXd &g_scale,
                           int new_first_axis, double weight_mod) {
    auto op_F = std::make_unique<Op_SAFE>(pdata, stim_thresh, weight_mod);
    op_F->safe_params.set_params(const_cast<double *>(tau1.data()), const_cast<double *>(tau2.data()),
                                 const_cast<double *>(tau3.data()), const_cast<double *>(a1.data()),
                                 const_cast<double *>(a2.data()), const_cast<double *>(a3.data()),
                                 const_cast<double *>(stim_limit.data()), const_cast<double *>(g_scale.data()));
    op_F->safe_params.swap_first_axes(new_first_axis);
    all_op.push_back(std::move(op_F));
}

void GroptParams::add_SAFE_vec(const Eigen::VectorXd &stim_thresh_vec, int new_first_axis, double weight_mod) {
    auto op_F = std::make_unique<Op_SAFE>(pdata, static_cast<int>(stim_thresh_vec.size()),
                                          const_cast<double *>(stim_thresh_vec.data()), weight_mod);
    op_F->safe_params.set_demo_params();
    op_F->safe_params.swap_first_axes(new_first_axis);
    all_op.push_back(std::move(op_F));
}

void GroptParams::add_SAFE_vec(const Eigen::VectorXd &stim_thresh_vec, const Eigen::VectorXd &tau1,
                               const Eigen::VectorXd &tau2, const Eigen::VectorXd &tau3, const Eigen::VectorXd &a1,
                               const Eigen::VectorXd &a2, const Eigen::VectorXd &a3, const Eigen::VectorXd &stim_limit,
                               const Eigen::VectorXd &g_scale, int new_first_axis, double weight_mod) {
    auto op_F = std::make_unique<Op_SAFE>(pdata, static_cast<int>(stim_thresh_vec.size()),
                                          const_cast<double *>(stim_thresh_vec.data()), weight_mod);
    op_F->safe_params.set_params(const_cast<double *>(tau1.data()), const_cast<double *>(tau2.data()),
                                 const_cast<double *>(tau3.data()), const_cast<double *>(a1.data()),
                                 const_cast<double *>(a2.data()), const_cast<double *>(a3.data()),
                                 const_cast<double *>(stim_limit.data()), const_cast<double *>(g_scale.data()));
    op_F->safe_params.swap_first_axes(new_first_axis);
    all_op.push_back(std::move(op_F));
}

void GroptParams::add_bvalue(double target, double tol, int start_idx0, int stop_idx0, double weight_mod, int mode,
                             double max_scale) {

    all_op.push_back(std::make_unique<Op_BValue>(pdata, target, tol, start_idx0, stop_idx0, weight_mod,
                                                 static_cast<BVALUE_MODE>(mode), max_scale));
}

void GroptParams::add_TV(double tv_lam, double weight_mod) {
    all_op.push_back(std::make_unique<Op_TV>(pdata, tv_lam, weight_mod));
}

void GroptParams::add_obj_identity(double weight_mod) {
    all_obj.push_back(std::make_unique<Op_Identity>(pdata, weight_mod));
}

void GroptParams::reset_op_weights() {
    for (auto &op : all_op) {
        // op->weight_mod = 1.0;
        op->spec_norm = 1.0;
        op->spec_norm2 = 1.0;
    }
    for (auto &op : all_obj) {
        // op->weight_mod = 1.0;
        op->spec_norm = 1.0;
        op->spec_norm2 = 1.0;
    }
}

double GroptParams::get_output_bvalue(const Eigen::VectorXd &X) {
    Op_BValue opB(pdata, 1, 1, -1, -1, 1, static_cast<BVALUE_MODE>(1), 1.01);
    opB.init();
    Eigen::VectorXd X_copy = X;
    return opB.get_bvalue(X_copy);
}

Eigen::VectorXd linear_interpolate(const Eigen::VectorXd &in, int out_size) {
    int in_size = in.size();
    if (out_size >= in_size) {
        return in;
    }

    Eigen::VectorXd out(out_size);
    double scale = static_cast<double>(in_size - 1) / (out_size - 1);

    for (int i = 0; i < out_size; ++i) {
        double in_idx_float = i * scale;
        int idx0 = static_cast<int>(floor(in_idx_float));
        int idx1 = idx0 + 1;

        if (idx1 >= in_size) { // Should only happen for the last element
            out(i) = in(in_size - 1);
        } else {
            double val0 = in(idx0);
            double val1 = in(idx1);
            double frac = in_idx_float - idx0;
            out(i) = val0 * (1.0 - frac) + val1 * frac;
        }
    }
    return out;
}

} // namespace Gropt

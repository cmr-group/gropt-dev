#include "spdlog/spdlog.h"

#include <map>

#include "gropt_params.hpp"

#include "op_bvalue.hpp"
#include "op_concomitant.hpp"
#include "op_diffbasin.hpp"
#include "op_eddy.hpp"
#include "op_gradient.hpp"
#include "op_identity.hpp"
#include "op_moment.hpp"
#include "op_safe.hpp"
#include "op_slew.hpp"
#include "op_tv.hpp"

namespace Gropt {

GroptParams::GroptParams() {}

void GroptParams::rebuild_fixer_from_set_vals() {
    pdata.fixer.setOnes(pdata.set_vals.size());
    for (int i = 0; i < pdata.set_vals.size(); i++) {
        if (!std::isnan(pdata.set_vals(i))) {
            pdata.fixer(i) = 0.0;
        }
    }
}

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

    rebuild_fixer_from_set_vals();

    pdata.X0.setOnes(N * Naxis);
    pdata.X0 *= .01;
    pdata.X0(0) = first_val;
    pdata.X0(N - 1) = last_val;

    vec_init_status = N;
}

void GroptParams::setvec_X0(const Eigen::VectorXd &_X0, int _Naxis, bool set_others) {
    int _N = static_cast<int>(_X0.size()) / _Naxis;

    if (_N > 0) {
        N = _N;
    }

    if (_Naxis > 0) {
        Naxis = _Naxis;
    }

    Ntot = N * Naxis;

    pdata.X0 = _X0;

    if (set_others) {
        pdata.inv_vec.setOnes(N * Naxis);

        pdata.set_vals.setZero(N * Naxis);
        pdata.set_vals.array() *= NAN;

        for (int j = 0; j < Naxis; j++) {
            pdata.set_vals((j * N)) = pdata.X0((j * N));
            pdata.set_vals((j * N) + N - 1) = pdata.X0((j * N) + N - 1);
        }

        rebuild_fixer_from_set_vals();
    }

    vec_init_status = N;
}

void GroptParams::setvec_set_vals(const Eigen::VectorXd &_set_vals, int _Naxis) {
    if (_Naxis <= 0) {
        throw std::invalid_argument("setvec_set_vals: Naxis must be > 0");
    }
    if (_set_vals.size() % _Naxis != 0) {
        throw std::invalid_argument("setvec_set_vals: set_vals.size() must be divisible by Naxis");
    }

    int _N = static_cast<int>(_set_vals.size()) / _Naxis;
    N = _N;
    Naxis = _Naxis;
    Ntot = N * Naxis;

    pdata.set_vals = _set_vals;
    rebuild_fixer_from_set_vals();

    if (pdata.X0.size() != Ntot) {
        pdata.X0.setOnes(Ntot);
        pdata.X0 *= .01;
    }
    for (int i = 0; i < pdata.set_vals.size(); i++) {
        if (!std::isnan(pdata.set_vals(i))) {
            pdata.X0(i) = pdata.set_vals(i);
        }
    }

    vec_init_status = N;
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

    spdlog::trace("diff_init N: {}", N);

    int ind_inv = (int)(TE / 2.0 / dt);
    pdata.inv_vec.setOnes(N);
    for (int i = ind_inv; i < N; i++) {
        pdata.inv_vec(i) = -1.0;
    }

    spdlog::trace("diff_init ind_inv: {}", ind_inv);

    int ind_90_end, ind_180_start, ind_180_end;
    ind_90_end = ceil(T_90 / dt);
    ind_180_start = floor((TE / 2.0 - T_180 / 2.0) / dt);
    ind_180_end = ceil((TE / 2.0 + T_180 / 2.0) / dt);

    spdlog::trace("diff_init inds: {}   {}   {}", ind_90_end, ind_180_start, ind_180_end);

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

    rebuild_fixer_from_set_vals();

    pdata.X0.setOnes(N);
    for (int i = 0; i < pdata.set_vals.size(); i++) {
        if (!isnan(pdata.set_vals(i))) {
            pdata.X0(i) = pdata.set_vals(i);
        } else {
            pdata.X0(i) = 1e-2; // Initial value for non-fixed points
        }
    }
    pdata.X0.array() *= pdata.inv_vec.array();

    vec_init_status = N;
}

int GroptParams::diff_init_preencode(double _dt, double _TE, double _T_90, double _T_180, double _T_readout,
                                     double _T_pre) {
    dt = _dt;
    Naxis = 1;

    double T_90 = _T_90;
    double T_180 = _T_180;
    double T_readout = _T_readout;
    double TE = _TE;

    int N_pre = floor(_T_pre / dt);

    N = N_pre + (int)((TE - T_readout) / dt) + 1;
    Ntot = N * Naxis;

    int ind_inv = N_pre + (int)(TE / 2.0 / dt);
    pdata.inv_vec.setOnes(N);
    for (int i = ind_inv; i < N; i++) {
        pdata.inv_vec(i) = -1.0;
    }

    int ind_90_end, ind_180_start, ind_180_end;
    ind_90_end = N_pre + ceil(T_90 / dt);
    ind_180_start = N_pre + floor((TE / 2.0 - T_180 / 2.0) / dt);
    ind_180_end = N_pre + ceil((TE / 2.0 + T_180 / 2.0) / dt);

    pdata.set_vals.setOnes(N);
    pdata.set_vals.array() *= NAN;
    for (int i = N_pre; i <= ind_90_end; i++) {
        pdata.set_vals(i) = 0.0;
    }
    for (int i = ind_180_start; i <= ind_180_end; i++) {
        pdata.set_vals(i) = 0.0;
    }
    pdata.set_vals(0) = 0.0;
    pdata.set_vals(N - 1) = 0.0;

    rebuild_fixer_from_set_vals();

    pdata.X0.setOnes(N);
    for (int i = 0; i < pdata.set_vals.size(); i++) {
        if (!isnan(pdata.set_vals(i))) {
            pdata.X0(i) = pdata.set_vals(i);
        } else {
            pdata.X0(i) = 1e-2; // Initial value for non-fixed points
        }
    }

    pdata.X0.array() *= pdata.inv_vec.array();

    vec_init_status = N;

    return N_pre;
}

void GroptParams::diff_init_deadtime(double _dt, double _TE, double _T_90, double _T_180, double _T_readout) {
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
    ind_180_end = ceil((TE / 2.0 + T_180 / 2.0) / dt);

    int live_time = N - ind_180_end;
    ind_180_start = ind_90_end + live_time - 1;

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

    rebuild_fixer_from_set_vals();

    pdata.X0.setOnes(N);
    for (int i = 0; i < pdata.set_vals.size(); i++) {
        if (!isnan(pdata.set_vals(i))) {
            pdata.X0(i) = pdata.set_vals(i);
        } else {
            pdata.X0(i) = 1e-2; // Initial value for non-fixed points
        }
    }

    pdata.X0.array() *= pdata.inv_vec.array();

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

    pdata.set_vals = set_vals_new;
    rebuild_fixer_from_set_vals();

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

    // Assign each operator a stable unique_name "<name>#<occurrence>" so duplicate-named operators
    // stay distinguishable (e.g. Moment#0, Moment#1). Used to match operators across solves for
    // warm starting (see warmstart.hpp); rebuilding the problem in the same order reproduces them.
    {
        std::map<std::string, int> name_counts;
        for (auto &op : all_op) {
            op->unique_name = op->name + "#" + std::to_string(name_counts[op->name]++);
        }
    }

    // Build the exact equality projector from any operators flagged use_projection (moments with
    // project=true, etc.). Each contributes its linear-equality
    // rows via append_eq_rows. If any participating operator's rows depend on the iterate
    // (eq_rows_vary -> relinearized), the solver rebuilds the projector each outer iteration.
    {
        bool any_proj = false;  // Are any projections on
        eq_proj_dynamic = false;  // Do the rows change
        for (auto &op : all_op) {
            if (op->use_projection) {
                any_proj = true;
                if (op->eq_rows_vary()) eq_proj_dynamic = true;
            }
        }
        if (any_proj) {
            build_eq_proj(pdata.X0, true); // initial build (logs rank/conditioning)
        } else {
            eq_proj.active = false;
        }
    }

    op_prep_status = N;

    spdlog::trace("GroptParams::prepare() end");
}

void GroptParams::add_gmax(double gmax, bool rot_variant, double weight_mod) {
    all_op.push_back(std::make_unique<Op_Gradient>(pdata, gmax, rot_variant, weight_mod));
}

void GroptParams::add_gmax_vec(const Eigen::VectorXd &gmax_vec, bool rot_variant, double weight_mod) {
    all_op.push_back(std::make_unique<Op_Gradient>(pdata, gmax_vec, rot_variant, weight_mod));
}

void GroptParams::add_smax(double smax, bool rot_variant, double weight_mod) {
    all_op.push_back(std::make_unique<Op_Slew>(pdata, smax, rot_variant, weight_mod));
}

void GroptParams::add_smax_vec(const Eigen::VectorXd &smax_vec, bool rot_variant, double weight_mod) {
    all_op.push_back(std::make_unique<Op_Slew>(pdata, smax_vec, rot_variant, weight_mod));
}

void GroptParams::build_eq_proj(const Eigen::VectorXd &x0, bool do_log) {
    // Collect linear-equality rows from every operator (linearized at x0), stack, and build.
    std::vector<Eigen::VectorXd> rows;
    std::vector<double> targets;
    for (auto &op : all_op) {
        op->append_eq_rows(rows, targets, x0);
    }

    int k = static_cast<int>(rows.size());
    if (k == 0) {
        eq_proj.active = false; // nothing to project this iteration (e.g. degenerate linearization)
        return;
    }

    Eigen::MatrixXd M(k, Ntot);
    Eigen::VectorXd t(k);
    for (int r = 0; r < k; r++) {
        M.row(r) = rows[r].transpose();
        t(r) = targets[r];
    }
    eq_proj.build(M, t, pdata.fixer, eq_proj_solver, eq_proj_rcond);

    if (do_log) {
        if (!eq_proj.rank_ok) {
            // TODO: This is wrong to some degree, it often works just fine when this ewrror is triggered.
            spdlog::warn("Equality projection: the free (non-fixed) DOFs cannot independently control the "
                          "{} projected row(s) -- nearly linearly dependent (scale-invariant cond = {:.2e}). "
                          "Free more of the waveform or reduce projected constraints.",
                          k, eq_proj.cond);
        } else {
            spdlog::info("Equality projection active: {} row(s), independence cond = {:.2e}", k, eq_proj.cond);
        }
    }
}

void GroptParams::add_concomitant(int start_idx, bool rot_variant, double weight_mod, double tol0,
                                  double target, bool fix_gamma, double gamma_fix, bool project,
                                  bool as_objective) {
    auto op = std::make_unique<Op_Concomitant>(pdata, start_idx, rot_variant, weight_mod, tol0, target);
    op->fix_gamma = fix_gamma;
    op->gamma_fix = gamma_fix;
    if (as_objective) {  // TODO: probably get rid of this, it didn't really help anything, though may be a good demo pathway for other operators
        // Augmented-Lagrangian objective: drives c(x)=0 from the all_obj path (CG LHS/RHS + scalar
        // dual). Not the optimization target, so it must not corrupt best-feasible scoring.
        op->is_score_obj = false;
        all_obj.push_back(std::move(op));
    } else {
        op->use_projection = project; // project=true: exact relinearized equality projection; else
        all_op.push_back(std::move(op)); // soft ADMM box/band constraint (prox)
    }
}

void GroptParams::add_moment(double order, double target, double tol0, std::string units, int moment_axis,
                             int start_idx0, int stop_idx0, int ref_idx0, double weight_mod, bool project,
                             bool absolute_tol) {
    auto op = std::make_unique<Op_Moment>(pdata, order, target, tol0, units, moment_axis, start_idx0, stop_idx0,
                                          ref_idx0, weight_mod);
    op->use_projection = project;      // exact null-space projection instead of ADMM penalty
    op->absolute_tol = absolute_tol;   // tol is an absolute per-order bound vs the default M0-anchored scaling
    all_op.push_back(std::move(op));
}

void GroptParams::add_SAFE(double stim_thresh, int new_first_axis, double weight_mod) {
    auto op_F = std::make_unique<Op_SAFE>(pdata, stim_thresh, weight_mod);
    op_F->safe_params.set_demo_params();
    op_F->safe_params.swap_first_axes(new_first_axis);
    op_F->safe_eps = safe_eps;
    all_op.push_back(std::move(op_F));
}

void GroptParams::add_SAFE(double stim_thresh, const Eigen::VectorXd &tau1, const Eigen::VectorXd &tau2,
                           const Eigen::VectorXd &tau3, const Eigen::VectorXd &a1, const Eigen::VectorXd &a2,
                           const Eigen::VectorXd &a3, const Eigen::VectorXd &stim_limit, const Eigen::VectorXd &g_scale,
                           int new_first_axis, double weight_mod) {
    auto op_F = std::make_unique<Op_SAFE>(pdata, stim_thresh, weight_mod);
    op_F->safe_params.set_params(tau1, tau2, tau3, a1, a2, a3, stim_limit, g_scale);
    op_F->safe_params.swap_first_axes(new_first_axis);
    op_F->safe_eps = safe_eps;
    all_op.push_back(std::move(op_F));
}

void GroptParams::add_SAFE_vec(const Eigen::VectorXd &stim_thresh_vec, int new_first_axis, double weight_mod) {
    auto op_F = std::make_unique<Op_SAFE>(pdata, stim_thresh_vec, weight_mod);
    op_F->safe_params.set_demo_params();
    op_F->safe_params.swap_first_axes(new_first_axis);
    op_F->safe_eps = safe_eps;
    all_op.push_back(std::move(op_F));
}

void GroptParams::add_SAFE_vec(const Eigen::VectorXd &stim_thresh_vec, const Eigen::VectorXd &tau1,
                               const Eigen::VectorXd &tau2, const Eigen::VectorXd &tau3, const Eigen::VectorXd &a1,
                               const Eigen::VectorXd &a2, const Eigen::VectorXd &a3, const Eigen::VectorXd &stim_limit,
                               const Eigen::VectorXd &g_scale, int new_first_axis, double weight_mod) {
    auto op_F = std::make_unique<Op_SAFE>(pdata, stim_thresh_vec, weight_mod);
    op_F->safe_params.set_params(tau1, tau2, tau3, a1, a2, a3, stim_limit, g_scale);
    op_F->safe_params.swap_first_axes(new_first_axis);
    op_F->safe_eps = safe_eps;
    all_op.push_back(std::move(op_F));
}

void GroptParams::add_bvalue(double target, double tol, int start_idx0, int stop_idx0, double weight_mod, int mode,
                             double max_scale, bool as_objective, bool linearize) {

    auto op = std::make_unique<Op_BValue>(pdata, target, tol, start_idx0, stop_idx0, weight_mod,
                                          static_cast<BVALUE_MODE>(mode), max_scale);
    if (as_objective) {
        op->linearize_obj = linearize;    // concave maximization -> linearize to RHS (DCA) by default
        all_obj.push_back(std::move(op)); // maximize b-value via the objective path (obj_weight = -1)
    } else {
        all_op.push_back(std::move(op)); // enforce b-value as a constraint
    }
}

void GroptParams::add_eddy(const Eigen::VectorXd &lam, double tol, double weight_mod, bool project) {
    auto op = std::make_unique<Op_Eddy>(pdata, lam, tol, weight_mod);
    op->use_projection = project; // exact null-space projection (eddy current == 0) instead of box ADMM
    all_op.push_back(std::move(op));
}

void GroptParams::add_TV(double tv_lam, double weight_mod, int order) {
    all_op.push_back(std::make_unique<Op_TV>(pdata, tv_lam, weight_mod, order));
}

void GroptParams::add_diff_basin(double window_time, double eps_factor, double gmax, double weight_mod,
                                 bool same_sign) {
    all_op.push_back(
        std::make_unique<Op_DiffBasin>(pdata, window_time, eps_factor, gmax, weight_mod, same_sign));
}

void GroptParams::add_obj_identity(double weight_mod) {
    all_obj.push_back(std::make_unique<Op_Identity>(pdata, weight_mod));
}

void GroptParams::reset_op_weights() {
    for (auto &op : all_op) {
        op->spec_norm = 1.0;
        op->spec_norm2 = 1.0;
    }
    for (auto &op : all_obj) {
        op->spec_norm = 1.0;
        op->spec_norm2 = 1.0;
    }
}

void GroptParams::print_op_details() {
    for (auto &op : all_op) {
        op->print_details();
    }
    for (auto &op : all_obj) {
        op->print_details();
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

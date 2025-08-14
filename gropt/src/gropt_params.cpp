#include "spdlog/spdlog.h"

#include "gropt_params.hpp"

#include "op_gradient.hpp"
#include "op_moment.hpp"
#include "op_slew.hpp"
#include "op_identity.hpp"
#include "op_bvalue.hpp"
#include "op_safe.hpp"
#include "op_tv.hpp"
#include "solver_groptsdmm.hpp"

namespace Gropt {

GroptParams::GroptParams() {
}

void GroptParams::vec_init_simple(int _N, int _Naxis, double first_val, double last_val) {
    if (_N > 0) {
        N = _N;
    }

    if (_Naxis > 0) {
        Naxis = _Naxis;
    }   

    Ntot = N * Naxis;

    inv_vec.setOnes(N * Naxis);
    
    set_vals.setZero(N * Naxis);
    set_vals.array() *= NAN;
    set_vals(0) = first_val;
    set_vals(N-1) = last_val;

    fixer.setOnes(N * Naxis);
    fixer(0) = 0.0;
    fixer(N-1) = 0.0;

    X0.setOnes(N * Naxis);
    X0 *= .01;
    X0(0) = first_val;
    X0(N-1) = last_val;

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

    X0.setOnes(N * Naxis);
    for (int i=0; i<N * Naxis; i++) {
        X0(i) = _X0[i];
    }

    if (set_others) {
        inv_vec.setOnes(N * Naxis);
        
        set_vals.setZero(N * Naxis);
        set_vals.array() *= NAN;

        fixer.setOnes(N * Naxis);

        for (int j = 0; j < Naxis; j++) {
            set_vals((j*N)) = X0((j*N));
            set_vals((j*N) + N-1) = X0((j*N) + N-1);

            fixer((j*N)) = 0.0;
            fixer((j*N) + N-1) = 0.0;
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

    N = (int)((TE-T_readout)/dt) + 1;
    Ntot = N * Naxis;

    int ind_inv = (int)(TE/2.0/dt);
    inv_vec.setOnes(N);
    for(int i = ind_inv; i < N; i++) {
        inv_vec(i) = -1.0;
    }

    int ind_90_end, ind_180_start, ind_180_end;
    ind_90_end = ceil(T_90/dt);
    ind_180_start = floor((TE/2.0 - T_180/2.0)/dt);
    ind_180_end = ceil((TE/2.0 + T_180/2.0)/dt);

    set_vals.setOnes(N);
    set_vals.array() *= NAN;
    for(int i = 0; i <= ind_90_end; i++) {
        set_vals(i) = 0.0;
    }
    for(int i = ind_180_start; i <= ind_180_end; i++) {
        set_vals(i) = 0.0;
    }
    set_vals(0) = 0.0;
    set_vals(N-1) = 0.0;


    fixer.setOnes(N);
    for(int i = 0; i < set_vals.size(); i++) {
        if (!isnan(set_vals(i))) {
            fixer(i) = 0.0;
        }
    }


    X0.setOnes(N);
    for(int i = 0; i < set_vals.size(); i++) {
        if (!isnan(set_vals(i))) {
            X0(i) = set_vals(i);
        } else {
            X0(i) = 1e-2;  // Initial value for non-fixed points
        }
    }

    vec_init_status = N;

}

void GroptParams::set_ils_solver(std::string _ils_method) 
{
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

    inv_vec.setOnes(N * Naxis);
    
    Eigen::VectorXd set_vals_new;
    set_vals_new.setZero(N * Naxis);
    set_vals_new.array() *= NAN;
    set_vals_new(0) = set_vals(0);
    set_vals_new(N-1) = set_vals(set_vals.size()-1);

    Eigen::VectorXd fixer_new;
    fixer_new.setOnes(N * Naxis);
    fixer_new(0) = 0.0;
    fixer_new(N-1) = 0.0;

    set_vals = set_vals_new;
    fixer = fixer_new;

    vec_init_status = N;
}

// This is a warm starter assuming that N has not changed
void GroptParams::warm_start_prev() {
    X0 = final_X;

    for (int i = 0; i < all_op.size(); i++) {
        all_op[i]->do_init_weights = false;
        all_op[i]->init();
        all_op[i]->reinit_parsdmm();
        all_op[i]->prep_parsdmm(X0);
    }
    for (int i = 0; i < all_obj.size(); i++) {
        all_obj[i]->init();
    }
}

void GroptParams::prepare() {

    spdlog::trace("GroptParams::init() start");

    if (N*Naxis != Ntot) {
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
        all_op[i]->init_parsdmm();
        all_op[i]->prep_parsdmm(X0);
    }
    for (int i = 0; i < all_obj.size(); i++) {
        all_obj[i]->init();
    }

    op_prep_status = N;

    spdlog::trace("GroptParams::init() end");
}

void GroptParams::add_gmax(double gmax, bool rot_variant, double weight_mod) {
    all_op.push_back(new Op_Gradient(*this, gmax, rot_variant, weight_mod));
}

void GroptParams::add_smax(double smax, bool rot_variant, double weight_mod) {
    all_op.push_back(new Op_Slew(*this, smax, rot_variant, weight_mod));
}

void GroptParams::add_moment(double order, double target, 
                             double tol0, std::string units, int moment_axis, 
                             int start_idx0, int stop_idx0, int ref_idx0, double weight_mod) {
    all_op.push_back(new Op_Moment(*this, order, target, tol0, units, 
                                   moment_axis, start_idx0, stop_idx0, ref_idx0, weight_mod));
}

void GroptParams::add_SAFE(double stim_thresh, 
                           double *tau1, double *tau2, double *tau3, 
                           double *a1, double *a2, double *a3,
                           double *stim_limit, double *g_scale,
                           int new_first_axis, bool demo_params, double weight_mod) 
{
    Op_SAFE* op_F = new Op_SAFE(*this, stim_thresh, weight_mod, false);
    if (demo_params) {
        op_F->safe_params.set_demo_params();
    } else {
        op_F->safe_params.set_params(tau1, tau2, tau3, a1, a2, a3, stim_limit, g_scale);
    }
    op_F->safe_params.swap_first_axes(new_first_axis);
    all_op.push_back(op_F);
}

void GroptParams::add_SAFE_vec(int N_vec, double *stim_thresh_vec,
                      double *tau1, double *tau2, double *tau3, 
                      double *a1, double *a2, double *a3,
                      double *stim_limit, double *g_scale,
                      int new_first_axis, bool demo_params, 
                      double weight_mod)
{
    Op_SAFE* op_F = new Op_SAFE(*this, N_vec, stim_thresh_vec, weight_mod, false);
    if (demo_params) {
        op_F->safe_params.set_demo_params();
    } else {
        op_F->safe_params.set_params(tau1, tau2, tau3, a1, a2, a3, stim_limit, g_scale);
    }
    op_F->safe_params.swap_first_axes(new_first_axis);
    all_op.push_back(op_F);
}

void GroptParams::add_bvalue(double target, double tol, int start_idx0, int stop_idx0, double weight_mod) {
    all_op.push_back(new Op_BValue(*this, target, tol, start_idx0, stop_idx0, weight_mod));
}

void GroptParams::add_TV(double tv_lam, double weight_mod)
{
    all_op.push_back(new Op_TV(*this, tv_lam, weight_mod));  
}

void GroptParams::add_obj_identity(double weight_mod) {
    all_obj.push_back(new Op_Identity(*this, weight_mod));
}

void GroptParams::solve() {
    if (op_prep_status != N) {
        spdlog::info("Operators do not seem prepared, calling prepare()");
        prepare();
    }

    SolverGroptSDMM solver;
    solver.solve(*this);
}


void GroptParams::test_reduce_and_solve() {
    SolverGroptSDMM solver;
    test_reduce_and_solve(solver);
}

void GroptParams::test_reduce_and_solve(SolverGroptSDMM solver) {
    vec_reduce_simple(1);
    X0 = linear_interpolate(final_X, N);


    for (int i = 0; i < all_op.size(); i++) {
        all_op[i]->init();
        all_op[i]->init_parsdmm();
        all_op[i]->prep_parsdmm(X0);
    }
    for (int i = 0; i < all_obj.size(); i++) {
        all_obj[i]->init();
    }
    op_prep_status = N;

    solver.solve(*this);
}



void GroptParams::solve(int min_iter, 
                        int max_iter, 
                        double gamma_x, 
                        double ils_tol, 
                        int ils_max_iter, 
                        int ils_min_iter, 
                        double ils_sigma,
                        double ils_tik_lam
                        ) 
{
    if (op_prep_status != N) {
        spdlog::info("Operators do not seem prepared, calling prepare()");
        prepare();
    }

    SolverGroptSDMM solver;
    solver.set_general_params(min_iter, max_iter, 20, gamma_x, 12000);
    solver.set_ils_params(ils_tol, ils_max_iter, ils_min_iter,
                          ils_sigma, ils_tik_lam);
    solver.solve(*this);
}

void GroptParams::get_output(double **out, int &out_size)
{
    out_size = final_X.size();
    *out = new double[out_size];
    for (int i = 0; i < out_size; i++) {
        (*out)[i] = final_X(i);
    }
}

Eigen::VectorXd linear_interpolate(const Eigen::VectorXd& in, int out_size) {
    int in_size = in.size();
    if (out_size >= in_size) {
        // Or throw an error, depending on desired behavior
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
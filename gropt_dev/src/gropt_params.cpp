#include "spdlog/spdlog.h"

#include "gropt_params.hpp"

#include "op_gradient.hpp"
#include "op_moment.hpp"
#include "op_slew.hpp"
#include "op_identity.hpp"
#include "op_bvalue.hpp"
#include "op_safe.hpp"
#include "solver_groptsdmm.hpp"

namespace Gropt {

GroptParams::GroptParams() {
}

void GroptParams::vec_init_simple() {
    inv_vec.setOnes(N * Naxis);
    
    set_vals.setZero(N * Naxis);
    set_vals.array() *= NAN;
    set_vals(0) = 0.0;
    set_vals(N-1) = 0.0;
    
    fixer.setOnes(N * Naxis);
    fixer(0) = 0.0;
    fixer(N-1) = 0.0;

    X0.setOnes(N * Naxis);
    X0 *= .01;
    X0(0) = 0.0;
    X0(N-1) = 0.0;

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

// This is a warm starter assuming that N ahs not changed
void GroptParams::warm_start_prev() {
    X0 = final_X;

    for (int i = 0; i < all_op.size(); i++) {
        all_op[i]->save_weights = true;
        all_op[i]->init();
        all_op[i]->reinit_parsdmm();
        all_op[i]->prep_parsdmm(X0);
    }
    for (int i = 0; i < all_obj.size(); i++) {
        all_obj[i]->init();
    }

}

void GroptParams::init() {

    spdlog::trace("GroptParams::init() start");

    if (vec_init_status != N) {
        spdlog::warn("set_vals and inv_vec were not initialized, calling vec_init_simple()");
        vec_init_simple();
    }

    for (int i = 0; i < all_op.size(); i++) {
        all_op[i]->init();
        all_op[i]->init_parsdmm();
        all_op[i]->prep_parsdmm(X0);
    }
    for (int i = 0; i < all_obj.size(); i++) {
        all_obj[i]->init();
    }
    spdlog::trace("GroptParams::init() end");
}

void GroptParams::add_gmax(double gmax) {
    all_op.push_back(new Op_Gradient(*this, gmax));
}

void GroptParams::add_smax(double smax) {
    all_op.push_back(new Op_Slew(*this, smax));
}

void GroptParams::add_moment(double order, double target) {
    all_op.push_back(new Op_Moment(*this, order, target));
}

void GroptParams::add_SAFE(double stim_thresh, 
                           double *tau1, double *tau2, double *tau3, 
                           double *a1, double *a2, double *a3,
                           double *stim_limit, double *g_scale,
                           int new_first_axis, bool demo_params) 
{
    Op_SAFE* op_F = new Op_SAFE(*this, stim_thresh);
    if (demo_params) {
        op_F->safe_params.set_demo_params();
    } else {
        op_F->safe_params.set_params(tau1, tau2, tau3, a1, a2, a3, stim_limit, g_scale);
    }
    op_F->safe_params.swap_first_axes(new_first_axis);
    all_op.push_back(op_F);
}

void GroptParams::add_bvalue(double target, double tol) {
    Op_BValue* op_Bval = new Op_BValue(*this);
    op_Bval->target = target;
    op_Bval->tol0 = tol;
    all_op.push_back(op_Bval);
}

void GroptParams::add_obj_identity(double weight_mod) {
    all_obj.push_back(new Op_Identity(*this, weight_mod));
}

void GroptParams::solve() {
    SolverGroptSDMM solver(*this);
    solver.solve();
}

void GroptParams::solve(int min_iter, 
                        int n_iter, 
                        double gamma_x, 
                        double ils_tol, 
                        int ils_max_iter, 
                        int ils_min_iter, 
                        double ils_sigma,
                        double ils_tik_lam
                        ) 
{
    SolverGroptSDMM solver(*this);
    solver.min_iter = min_iter;
    solver.N_iter = n_iter;
    solver.gamma_x = gamma_x;
    solver.ils_tol = ils_tol;
    solver.ils_max_iter = ils_max_iter;
    solver.ils_min_iter = ils_min_iter;
    solver.ils_sigma = ils_sigma;
    solver.ils_tik_lam = ils_tik_lam;
    solver.solve();
}

void GroptParams::get_output(double **out, int &out_size)
{
    out_size = final_X.size();
    *out = new double[out_size];
    for (int i = 0; i < out_size; i++) {
        (*out)[i] = final_X(i);
    }
}

} // namespace Gropt
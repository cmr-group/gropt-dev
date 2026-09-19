#ifndef GROPT_PARAMS_H
#define GROPT_PARAMS_H

#include "Eigen/Dense"
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "op_main.hpp"
#include "problem_data.hpp"

namespace Gropt {

enum ILSMethod {
    CG,
    NLCG,
    BiCGstabl,
};

// Linear solver for the equality projection (project=true constraints):
//   EQ_LDLT: LDLT of the unit-row Gram Mhat Mhatᵀ; fast, but fails on singular or collinear rows.
//   EQ_COD:  complete-orthogonal decomposition of Mhat (condition number not squared); rank-revealing, so
//            near-collinear rows (e.g. several eddy time constants) still give an idempotent projector.
enum EqProjSolver {
    EQ_LDLT,
    EQ_COD,
};

class Operator;

struct SolveResult {
    Eigen::VectorXd X;      // solution waveform (SDMM: best feasible iterate, else the last one)
    bool converged = false; // every constraint feasible at X
    int n_iter = 0;         // outer iterations
    int n_feval = 0;        // total inner linear-solver iterations
    double dt = 0.0;
    double bvalue = 0.0;    // b-value of X if a b-value operator is present, else 0
};

// Exact projection onto {x : M x = t} for linear-equality constraints flagged project=true.
// Operators supply their rows via Operator::append_eq_rows.
//   M    = stacked equality rows (k x Ntot)
//   Mtil = M * diag(fixer), i.e. M with fixed DOFs zeroed
//   Mhat = D^-1 Mtil, rows scaled to unit norm for scale-invariant conditioning and rank checks
// Corrections lie in range(Mhatᵀ), so only free DOFs move; the residual M x - t uses the full M so
// nonzero fixed samples are accounted for.
struct EqualityProjection {
    bool active = false;
    EqProjSolver solver = EQ_LDLT; // which factorization the project_* methods use
    Eigen::MatrixXd M;        // k x Ntot, full equality rows (for the residual M x - t)
    Eigen::MatrixXd Mhat;     // k x Ntot, free-masked rows scaled to unit norm (D^-1 * Mtil)
    Eigen::VectorXd inv_rownorm; // k, 1/||free-masked row|| (the D^-1 scaling; 0 for fully-fixed rows)
    Eigen::VectorXd t;        // k, targets
    Eigen::LDLT<Eigen::MatrixXd> G;                              // EQ_LDLT: factor of the unit-row Gram
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod; // EQ_COD: rank-revealing factor of Mhat
    double cond = 0.0;        // condition number of the unit-row Gram (scale-invariant; -1 => singular)
    int n_uncontrolled = 0;   // rows with no free samples (the projection cannot enforce them)
    bool ill_conditioned = false; // Gram cond >= 1e15 or singular: the LDLT solve loses accuracy

    void build(const Eigen::MatrixXd &M_in, const Eigen::VectorXd &t_in, const Eigen::VectorXd &fixer,
               EqProjSolver solver_in, double cod_rcond) {
        solver = solver_in;
        M = M_in;
        t = t_in;
        Eigen::MatrixXd Mtil = M_in * fixer.asDiagonal();

        Eigen::VectorXd rownorm = Mtil.rowwise().norm();
        inv_rownorm.setZero(rownorm.size());
        Mhat = Mtil;
        for (int i = 0; i < Mhat.rows(); i++) {
            if (rownorm(i) > 0.0) {
                inv_rownorm(i) = 1.0 / rownorm(i);
                Mhat.row(i) *= inv_rownorm(i);
            } else {
                Mhat.row(i).setZero(); // fully-fixed row: free DOFs cannot control it
            }
        }

        // Conditioning of the unit-row Gram over rows that touch free samples (the rest are counted in
        // n_uncontrolled). Moment rows (t^0..t^k) are Hilbert-like, so cond grows quickly with order.
        Eigen::MatrixXd Ghat = Mhat * Mhat.transpose();
        n_uncontrolled = static_cast<int>((rownorm.array() == 0.0).count());
        std::vector<int> live;
        for (int i = 0; i < rownorm.size(); i++) {
            if (rownorm(i) > 0.0) live.push_back(i);
        }
        cond = 1.0;
        ill_conditioned = false;
        if (live.size() > 1) {
            Eigen::MatrixXd Glive(live.size(), live.size());
            for (size_t a = 0; a < live.size(); a++) {
                for (size_t b = 0; b < live.size(); b++) Glive(a, b) = Ghat(live[a], live[b]);
            }
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Glive);
            double lo = es.eigenvalues().minCoeff();
            double hi = es.eigenvalues().maxCoeff();
            cond = (lo > 0.0) ? (hi / lo) : -1.0; // -1 flags singular
            ill_conditioned = (lo <= 1e-15 * hi);
        }

        // Factorize (see EqProjSolver). cod_rcond is the EQ_COD rank threshold; <= 0 keeps Eigen's default.
        if (solver == EQ_COD) {
            cod.compute(Mhat);
            if (cod_rcond > 0.0) cod.setThreshold(cod_rcond);
        } else {
            G.compute(Ghat);
        }
        active = (M.rows() > 0);
    }

    // ILS_CG calls project_affine at the start of each inner solve and project_dir on the residual and
    // every A*p, so x stays feasible throughout CG.

    // Move x onto {M x = t}: subtract the min-norm free-DOF correction Δ with Mhat Δ = D^-1 (M x - t).
    void project_affine(Eigen::VectorXd &x) const {
        if (!active) return;
        Eigen::VectorXd r = (M * x - t).cwiseProduct(inv_rownorm); // D^-1 (M x - t)
        if (solver == EQ_COD) {
            x.noalias() -= cod.solve(r);
        } else {
            x.noalias() -= Mhat.transpose() * G.solve(r);
        }
    }

    // Project a direction onto null(Mhat), removing its free row-space component.
    void project_dir(Eigen::VectorXd &v) const {
        if (!active) return;
        if (solver == EQ_COD) {
            v.noalias() -= cod.solve(Mhat * v);
        } else {
            v.noalias() -= Mhat.transpose() * G.solve(Mhat * v);
        }
    }
};

class GroptParams {

  public:
    ProblemData pdata;

    // References into pdata
    double &dt = pdata.dt;
    int &N = pdata.N;
    int &Naxis = pdata.Naxis;

    int Ntot = 10;

    int vec_init_status = -1;
    int op_prep_status = -1;
    size_t op_prep_count = 0; // all_op + all_obj size at the last prepare()

    std::vector<std::unique_ptr<Operator>> all_op;
    std::vector<std::unique_ptr<Operator>> all_obj;

    // Scale linearized (DCA) objective pulls to magnitude |obj_weight|, independent of ||AᵀA x||.
    bool normalize_obj = false;

    // SAFE softabs smoothing [T/m/s] (see Op_SAFE::safe_eps); copied into each Op_SAFE by add_SAFE, so set it first.
    double safe_eps = 0.0;

    // Equality projector for project=true constraints, built in prepare(). eq_proj_dynamic: some projected
    // operator has iterate-dependent rows, so the solver rebuilds it every outer iteration.
    EqualityProjection eq_proj;
    bool eq_proj_dynamic = false;

    // Factorization used by eq_proj (see EqProjSolver).
    EqProjSolver eq_proj_solver = EQ_LDLT;
    // EQ_COD rank tolerance: singular values below eq_proj_rcond * largest are treated as zero (larger
    // drops more near-dependent rows); <= 0 uses Eigen's default. Ignored by EQ_LDLT.
    double eq_proj_rcond = 1e-10;

    // (Re)build eq_proj from every operator's append_eq_rows, linearized at x0; do_log logs
    // conditioning warnings.
    void build_eq_proj(const Eigen::VectorXd &x0, bool do_log = false);

    ILSMethod ils_method = CG;

    GroptParams();
    ~GroptParams() = default;

    // Move-only due to unique_ptr members
    GroptParams(GroptParams &&) = default;
    GroptParams &operator=(GroptParams &&) = default;
    GroptParams(const GroptParams &) = delete;
    GroptParams &operator=(const GroptParams &) = delete;

    void vec_init_simple(int _N, int _Naxis, double first_val, double last_val);
    void diff_init(double _dt, double _TE, double _T_90, double _T_180, double _T_readout);
    int diff_init_preencode(double _dt, double _TE, double _T_90, double _T_180, double _T_readout, double _T_pre);
    void diff_init_deadtime(double _dt, double _TE, double _T_90, double _T_180, double _T_readout);
    void setvec_X0(const Eigen::VectorXd &_X0, int _Naxis, bool set_others);
    void setvec_set_vals(const Eigen::VectorXd &_set_vals, int _Naxis);

    void vec_reduce_simple(int N_reduce);

    void prepare();
    // True if N changed or operators were added since the last prepare()
    bool needs_prepare() const {
        return op_prep_status != N || op_prep_count != all_op.size() + all_obj.size();
    }
    // Warm-start key ("<name>#<occurrence>") of each constraint operator, in all_op order. prepare() copies
    // these into Operator::unique_name; this call does not prepare.
    std::vector<std::string> get_op_keys() const;

    void set_ils_solver(std::string ils_method);

    void add_gmax(double gmax, bool rot_variant, double weight_mod);
    void add_gmax_vec(const Eigen::VectorXd &gmax_vec, bool rot_variant, double weight_mod);
    void add_smax(double smax, bool rot_variant, double weight_mod);
    void add_smax_vec(const Eigen::VectorXd &smax_vec, bool rot_variant, double weight_mod);
    // rot_variant is currently ignored: the energy balance is always computed across all axes.
    void add_concomitant(int start_idx, bool rot_variant, double weight_mod, double tol0 = 0.1,
                         double target = 1.0, bool project = false);
    void add_moment(double order, double target, double tol0, std::string units, int moment_axis, int start_idx0,
                    int stop_idx0, int ref_idx0, double weight_mod, bool project = false, bool absolute_tol = false);

    void add_SAFE(double stim_thresh, int new_first_axis, double weight_mod);
    void add_SAFE(double stim_thresh, const Eigen::VectorXd &tau1, const Eigen::VectorXd &tau2,
                  const Eigen::VectorXd &tau3, const Eigen::VectorXd &a1, const Eigen::VectorXd &a2,
                  const Eigen::VectorXd &a3, const Eigen::VectorXd &stim_limit, const Eigen::VectorXd &g_scale,
                  int new_first_axis, double weight_mod);

    void add_SAFE_vec(const Eigen::VectorXd &stim_thresh_vec, int new_first_axis, double weight_mod);
    void add_SAFE_vec(const Eigen::VectorXd &stim_thresh_vec, const Eigen::VectorXd &tau1, const Eigen::VectorXd &tau2,
                      const Eigen::VectorXd &tau3, const Eigen::VectorXd &a1, const Eigen::VectorXd &a2,
                      const Eigen::VectorXd &a3, const Eigen::VectorXd &stim_limit, const Eigen::VectorXd &g_scale,
                      int new_first_axis, double weight_mod);

    void add_eddy(const Eigen::VectorXd &lam, double tol, double weight_mod, bool project = false);

    void add_bvalue(double target, double tol, int start_idx0, int stop_idx0, double weight_mod, int mode,
                    double max_scale, bool as_objective = false, bool linearize = true);
    void add_TV(double tv_lam, double weight_mod, int order = 1);

    void add_diff_basin(double window_time, double eps_factor, double gmax, double weight_mod = 1.0,
                        bool same_sign = false);

    void add_obj_identity(double weight_mod);

    void reset_op_weights();

    void print_op_details();

    double get_output_bvalue(const Eigen::VectorXd &X);

  private:
    // pdata.fixer from set_vals: 1.0 where free (NaN), 0.0 where fixed.
    void rebuild_fixer_from_set_vals();
};

Eigen::VectorXd linear_interpolate(const Eigen::VectorXd &in, int out_size);

} // namespace Gropt

#endif

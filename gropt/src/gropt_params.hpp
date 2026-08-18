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
//   EQ_LDLT -- exact LDLT of the unit-row Gram. Fast and precise for well-conditioned constraint sets
//              (e.g. moments + concomitant). Blows up if rows are singular / collinear.
//   EQ_COD  -- rank-revealing complete-orthogonal decomposition of Mhat ITSELF (not the Gram, so the
//              condition number is not squared). Handles singular / near-collinear rows (e.g. multiple
//              eddy time-constants) and stays a TRUE idempotent projector, so it is safe inside the CG.
enum EqProjSolver {
    EQ_LDLT,
    EQ_COD,
};

class Operator; // Forward declaration of Operator class

struct SolveResult {
    Eigen::VectorXd X;
    bool converged = false;
    int n_iter = 0;
    int n_feval = 0;
    double dt = 0.0;
    double bvalue = 0.0;
};

// Exact null-space projection for linear-equality constraints flagged project=true.  
// Operators supply their rows via Operator::append_eq_rows.
//   M    = stack of equality functionals (k x Ntot)
//   Mtil (M~) = M with FIXED DOFs zeroed (free-mask = pdata.fixer): M * diag(fixer)
//   Mhat (M^) = row-scaled Mtil to unit norm (D^-1 * Mtil), for scale-invariant conditioning and rank-checking
//   Gram = Gram matrix
// The projection only moves FREE DOFs (Mtilᵀ is zero on fixed rows), and the residual (M x - t)
// uses the full M so nonzero fixed values' contributions are accounted for.
struct EqualityProjection {
    bool active = false;
    EqProjSolver solver = EQ_LDLT; // which factorization the project_* methods use
    Eigen::MatrixXd M;        // k x Ntot, full moment rows (for the affine residual M x - t)
    Eigen::MatrixXd Mhat;     // k x Ntot, free-masked rows scaled to UNIT norm (D^-1 * Mtil)
    Eigen::VectorXd inv_rownorm; // k, 1/||free-masked row|| (the D^-1 scaling; 0 for fully-fixed rows)
    Eigen::VectorXd t;        // k, targets
    Eigen::LDLT<Eigen::MatrixXd> G;                              // EQ_LDLT: exact factor of the unit-row Gram
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod; // EQ_COD: rank-revealing factor of Mhat
    double cond = 0.0;        // condition number of the unit-row Gram (scale-invariant; -1 => singular)
    bool rank_ok = true;      // false if the free DOFs can't independently control all rows

    void build(const Eigen::MatrixXd &M_in, const Eigen::VectorXd &t_in, const Eigen::VectorXd &fixer,
               EqProjSolver solver_in, double cod_rcond) {
        solver = solver_in;
        M = M_in;
        t = t_in;
        Eigen::MatrixXd Mtil = M_in * fixer.asDiagonal();

        // Row-scale to unit norm (Mhat).
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

        // Rank/conditioning check on the SCALE-INVARIANT unit-row Gram (matrix of cosine angles
        // between constraint rows). A near-zero eigenvalue here means two constraints are nearly
        // linearly dependent / the free DOFs genuinely can't separate them -- a real ill-posedness.
        Eigen::MatrixXd Ghat = Mhat * Mhat.transpose();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ghat);
        double lo = es.eigenvalues().minCoeff();
        double hi = es.eigenvalues().maxCoeff();
        cond = (lo > 0.0) ? (hi / lo) : -1.0; // -1 flags singular
        rank_ok = (rownorm.minCoeff() > 0.0) && (hi > 0.0) && (lo >= 1e-12 * hi);

        // Factorize per the chosen solver:
        //  EQ_LDLT: exact LDLT of the unit-row Gram -- fast and precise for well-conditioned constraint
        //           sets, but fails (blows up) if rows are singular / collinear.
        //  EQ_COD : Complete-orthogonal decomposition of Mhat ITSELF (not the Gram), so
        //           the condition number is not squared -> better precision on ill-conditioned moments.
        //           cod.solve gives the min-norm solution, which is a TRUE projector 
        //           cod_rcond is the rank tolerance: singular values below cod_rcond times
        //           the largest are treated as zero (<=0 uses Eigen's default).
        if (solver == EQ_COD) {
            cod.compute(Mhat);
            if (cod_rcond > 0.0) cod.setThreshold(cod_rcond);
        } else {
            G.compute(Ghat);
        }
        active = (M.rows() > 0);
    }

    // project_affine = move the point onto the surface
    //     used to re-initialize the starting vector, and before every CG iteration.
    // project_dir = remove any part of this arrow that points off the surface.
    //     used within the CG iterations to make sure directions always stay within the equality constraint manifold. 
    // A feasible starting point plus only-null-space updates means x stays feasible for free through all of CG


    // Force the TOTAL moment to t, moving only free DOFs. Min-norm correction Δ with Mhat Δ = D⁻¹(Mx-t):
    //   EQ_LDLT: Δ = Mhatᵀ (Mhat Mhatᵀ)⁻¹ r   EQ_COD: Δ = cod.solve(r) (min-norm, rank-robust)
    void project_affine(Eigen::VectorXd &x) const {
        if (!active) return;
        Eigen::VectorXd r = (M * x - t).cwiseProduct(inv_rownorm); // D^-1 (M x - t)
        if (solver == EQ_COD) {
            x.noalias() -= cod.solve(r);
        } else {
            x.noalias() -= Mhat.transpose() * G.solve(r);
        }
    }

    // Remove the (free) row-space component of a direction (project onto null(Mhat)).
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

    // Convenience accessors that forward to pdata
    double &dt = pdata.dt;
    int &N = pdata.N;
    int &Naxis = pdata.Naxis;

    int Ntot = 10;

    int vec_init_status = -1;
    int op_prep_status = -1;

    std::vector<std::unique_ptr<Operator>> all_op;
    std::vector<std::unique_ptr<Operator>> all_obj;

    // If true, self-normalize the objective so weight_mod sets a CONSTANT step magnitude (the pull no longer grows with ||AᵀA x||).
    bool normalize_obj = false;

    // Softabs smoothing (slew units) applied to every SAFE op's |.| (copied to each Op_SAFE in add_SAFE).
    // 0 = hard abs (original). >0 smooths the sign through zero -> removes the near-zero sign churn that
    // makes SAFE far more unstable than the linear eddy constraint. See Op_SAFE::safe_eps.
    double safe_eps = 0.0;

    // Exact null-space projection for linear-equality constraints flagged project=true. Built in
    // prepare(); if any participating operator has iterate-dependent rows (eq_proj_dynamic), the
    // solver rebuilds it each outer iteration via build_eq_proj(X).
    EqualityProjection eq_proj;
    bool eq_proj_dynamic = false;

    // Which linear solver the equality projection uses.
    EqProjSolver eq_proj_solver = EQ_LDLT;
    // Rank tolerance for EQ_COD only (singular values below eq_proj_rcond * largest are dropped);
    // ignored by EQ_LDLT. Larger = drop more near-dependent rows; <= 0 uses Eigen's default threshold.
    double eq_proj_rcond = 1e-10;

    // (Re)build eq_proj from every operator's append_eq_rows, linearized at x0. do_log emits the
    // rank/conditioning message.
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

    void set_ils_solver(std::string ils_method);

    void add_gmax(double gmax, bool rot_variant, double weight_mod);
    void add_gmax_vec(const Eigen::VectorXd &gmax_vec, bool rot_variant, double weight_mod);
    void add_smax(double smax, bool rot_variant, double weight_mod);
    void add_smax_vec(const Eigen::VectorXd &smax_vec, bool rot_variant, double weight_mod);
    void add_concomitant(int start_idx, bool rot_variant, double weight_mod, double tol0 = 0.1,
                         double target = 1.0, bool fix_gamma = false, double gamma_fix = 1.0,
                         bool project = false, bool as_objective = false);
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
    // Rebuilds pdata.fixer from the NaN pattern of pdata.set_vals:
    // 1.0 where free (NaN), 0.0 where fixed (finite).
    void rebuild_fixer_from_set_vals();
};

Eigen::VectorXd linear_interpolate(const Eigen::VectorXd &in, int out_size);

} // namespace Gropt

#endif

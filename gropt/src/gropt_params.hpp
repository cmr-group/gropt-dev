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

class Operator; // Forward declaration of Operator class

struct SolveResult {
    Eigen::VectorXd X;
    bool converged = false;
    int n_iter = 0;
    int n_feval = 0;
    double dt = 0.0;
    double bvalue = 0.0;
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
    void add_smax(double smax, bool rot_variant, double weight_mod);
    void add_smax_vec(const Eigen::VectorXd &smax_vec, bool rot_variant, double weight_mod);
    void add_concomitant(int start_idx, bool rot_variant, double weight_mod);
    void add_moment(double order, double target, double tol0, std::string units, int moment_axis, int start_idx0,
                    int stop_idx0, int ref_idx0, double weight_mod);

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

    void add_eddy(const Eigen::VectorXd &lam, double tol, double weight_mod);

    void add_bvalue(double target, double tol, int start_idx0, int stop_idx0, double weight_mod, int mode,
                    double max_scale);
    void add_TV(double tv_lam, double weight_mod);

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

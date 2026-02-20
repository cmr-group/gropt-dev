#ifndef GROPT_PARAMS_H
#define GROPT_PARAMS_H

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include "Eigen/Dense"

#include "problem_data.hpp"
#include "op_main.hpp"

namespace Gropt {

enum ILSMethod {
  CG,
  NLCG,
  BiCGstabl,
};

class Operator;  // Forward declaration of Operator class

struct SolveResult {
    Eigen::VectorXd X;
    bool converged = false;
    int n_iter = 0;
    int n_feval = 0;
};

class GroptParams
{

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
        GroptParams(GroptParams&&) = default;
        GroptParams& operator=(GroptParams&&) = default;
        GroptParams(const GroptParams&) = delete;
        GroptParams& operator=(const GroptParams&) = delete;

        void vec_init_simple(int _N, int _Naxis, double first_val, double last_val);
        void diff_init(double _dt, double _TE, double _T_90, double _T_180, double _T_readout);
        void setvec_X0(int _N, int _Naxis, double *_X0, bool set_others);

        void vec_reduce_simple(int N_reduce);

        void prepare();

        void set_ils_solver(std::string ils_method);

        void add_gmax(double gmax, bool rot_variant, double weight_mod);
        void add_smax(double smax, bool rot_variant, double weight_mod);
        void add_moment(double order, double target, double tol0, std::string units,
                        int moment_axis, int start_idx0, int stop_idx0, int ref_idx0, double weight_mod);
        void add_SAFE(double stim_thresh,
                      double *tau1, double *tau2, double *tau3,
                      double *a1, double *a2, double *a3,
                      double *stim_limit, double *g_scale,
                      int new_first_axis, bool demo_params,
                      double weight_mod);
        void add_SAFE_vec(int N_vec, double *stim_thresh_vec,
                      double *tau1, double *tau2, double *tau3,
                      double *a1, double *a2, double *a3,
                      double *stim_limit, double *g_scale,
                      int new_first_axis, bool demo_params,
                      double weight_mod);

        void add_bvalue(double target, double tol, int start_idx0, int stop_idx0,
                        double weight_mod, int mode, double max_scale);
        void add_TV(double tv_lam, double weight_mod);

        void add_obj_identity(double weight_mod);

        double get_output_bvalue(const Eigen::VectorXd &X);
};

Eigen::VectorXd linear_interpolate(const Eigen::VectorXd& in, int out_size);

} // namespace Gropt


#endif

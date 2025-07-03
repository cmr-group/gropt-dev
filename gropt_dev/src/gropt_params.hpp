#ifndef GROPT_PARAMS_H
#define GROPT_PARAMS_H

#include <iostream> 
#include <string>
#include <vector>
#include "Eigen/Dense"

#include "op_main.hpp"

namespace Gropt {

enum SolverMethod {
  GROPT_SDMM,
}; 

enum ILSMethod {
  CG,
}; 

class Operator;  // Forward declaration of Operator class
class GroptParams
{
    public:
        std::vector<Operator*> all_op;
        std::vector<Operator*> all_obj;

        Eigen::VectorXd X0;
        Eigen::VectorXd final_X;

        Eigen::VectorXd inv_vec;
        Eigen::VectorXd set_vals;
        Eigen::VectorXd fixer;

        ILSMethod ils_method = CG;
        SolverMethod solver_method = GROPT_SDMM;

        double dt = 10e-6;
        int N = 10;
        int Naxis = 1;

        int iiter;
        int final_good = 0;

        int vec_init_status = -1;

        GroptParams();
        ~GroptParams() = default;

        void vec_init_simple();
        void diff_init(double _dt, double _TE, double _T_90, double _T_180, double _T_readout);

        void init();
        void warm_start_prev();

        void add_gmax(double gmax);
        void add_smax(double smax);
        void add_moment(double order, double target);
        void add_SAFE(double stim_thresh);
        void add_bvalue(double target, double tol);

        void add_obj_identity(double weight_mod);

        void solve();
        void solve(int min_iter, 
                   int n_iter, 
                   double gamma_x, 
                   double ils_tol, 
                   int ils_max_iter, 
                   int ils_min_iter, 
                   double ils_sigma);

        void get_output(double **out, int &out_size);
};

} // namespace Gropt


#endif
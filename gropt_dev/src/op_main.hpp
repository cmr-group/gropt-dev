#ifndef OP_MAIN_H
#define OP_MAIN_H

/**
 * This is our main parent class for every constraint and regularization term in GrOpt
 * The main functions that we use here are the reweighting schemes, where most operators will
 * then just need to implement the forward, transpose, and proximal mapping operators.
 */

#include <iostream> 
#include <string>
#include <vector>
#include <math.h>  
#include "Eigen/Dense"

#include "gropt_params.hpp"

namespace Gropt {

class GroptParams;  // Forward declaration of GroptParams class
class Operator  // This is the main parent class for every operator in GrOpt 
{  
    public:   
        std::string name;

        GroptParams *gparams;
        
        int N;
        int Naxis;
        int Ntot;
        double dt;
        int Ax_size;

        bool rot_variant = true; 
        bool do_init_weights = true;
        
        double target;
        double tol0;
        double tol;
        double cushion = 1e-2;  // This is the cushion factor, which is used to reduce the tolerance by a factor of 1-cushion
        
        double spec_norm2;
        double spec_norm;

        double weight = 1.0;
        double obj_weight = 1.0;
        double weight_mod = 1.0;

        double gamma = 1.5;

        

        Eigen::VectorXd x_temp;
        Eigen::VectorXd x_temp_obj;
        Eigen::VectorXd Ax_temp;

        double feas_check; 
        double r_feas;
        Eigen::VectorXd feas_temp;

        std::vector<int> hist_feas;
        std::vector<double> hist_r_feas;

        // ----------------------------------------
        // PAR-SDMM Specific Variables
        // ---------------------------------------- 
        bool do_rw = true;
        bool do_gamma = true;
        bool do_weight = true;
        bool do_scalelim = true; 

        // Some of these should move to solver settings
        int rw_interval = 8;
        double e_corr = 0.4;
        double rw_eps = 1e-36;
        double rw_scalelim = 1.5;

        Eigen::VectorXd y0;
        Eigen::VectorXd y1;
        Eigen::VectorXd z0;
        Eigen::VectorXd z1;

        Eigen::VectorXd s0;
        Eigen::VectorXd s1;

        Eigen::VectorXd yhat1;
        Eigen::VectorXd dyhat;
        Eigen::VectorXd dy;
        Eigen::VectorXd dhhat;
        Eigen::VectorXd dghat;

        Eigen::VectorXd yhat00;
        Eigen::VectorXd y00;
        Eigen::VectorXd s00;
        Eigen::VectorXd z00;

        // ----------------------------------------
        // OSQP Specific Variables
        // ---------------------------------------- 
        // Eigen::VectorXd xhatk1;
        // Eigen::VectorXd xk0;
        // Eigen::VectorXd xk1;
        // Eigen::VectorXd r_dual; 

        // Eigen::VectorXd yk0;
        // Eigen::VectorXd yk1;
        // Eigen::VectorXd zk0;
        // Eigen::VectorXd zk1;
        // Eigen::VectorXd zhatk1;
        // Eigen::VectorXd dyk;
        // Eigen::VectorXd dzk;
        // Eigen::VectorXd r_primal;

        // std::vector<Eigen::VectorXd> hist_obj;
        // std::vector<Eigen::VectorXd> hist_Ax;
        // std::vector<Eigen::VectorXd> hist_y;
        // std::vector<Eigen::VectorXd> hist_z;
        // std::vector<Eigen::VectorXd> hist_Aty;


        // ----------------------------------------
        Operator() = default;
        Operator(GroptParams &_gparams);
        virtual ~Operator();

        virtual void init();
        virtual void init_parsdmm();
        virtual void reinit_parsdmm();

        virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void forward_op(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose_op(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void check(Eigen::VectorXd &X);
        virtual void prox(Eigen::VectorXd &X);
        virtual void get_feas(Eigen::VectorXd &s);
        virtual void prep_parsdmm(Eigen::VectorXd &X);
        virtual void reweight_parsdmm();
        virtual void add_Atb(Eigen::VectorXd &b);
        void add_AtAx(Eigen::VectorXd &x, Eigen::VectorXd &out);
        virtual void add_obj(Eigen::VectorXd &x, Eigen::VectorXd &out);
};

}  // close "namespace Gropt"

#endif
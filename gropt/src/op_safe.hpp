#ifndef OP_SAFE_H
#define OP_SAFE_H

/**
 * Constriant on FASE model prediction of waveforms.
 */

#include <iostream> 
#include <string>
#include <math.h>  
#include "Eigen/Dense"

#include "op_main.hpp"

namespace Gropt {

class SAFEParams {

    public:  // Code is too messy if this isnt public, any better way?

        std::vector<double> tau1 = std::vector<double>(3, 0.0);
        std::vector<double> tau2 = std::vector<double>(3, 0.0);
        std::vector<double> tau3 = std::vector<double>(3, 0.0);
        std::vector<double> a1 = std::vector<double>(3, 0.0);
        std::vector<double> a2 = std::vector<double>(3, 0.0);
        std::vector<double> a3 = std::vector<double>(3, 0.0);
        std::vector<double> stim_limit = std::vector<double>(3, 0.0);
        std::vector<double> g_scale = std::vector<double>(3, 0.0);

        // These aren't real parameters because they depend on dt, maybe they should be in Op_SAFE?
        std::vector<double> alpha1 = std::vector<double>(3, 0.0);
        std::vector<double> alpha2 = std::vector<double>(3, 0.0);
        std::vector<double> alpha3 = std::vector<double>(3, 0.0);

    
        SAFEParams() = default;
        void set_demo_params();
        void set_params(double *_tau1, double *_tau2, double *_tau3, 
                            double *_a1, double *_a2, double *_a3,
                            double *_stim_limit, double *_g_scale);
        void calc_alphas(double dt);
        void swap_first_axes(int new_first_axis);
};

class Op_SAFE : public Operator
{  
    protected:
        double stim_thresh;
        bool true_safe;

        Eigen::VectorXd signs;
        Eigen::VectorXd stim1;
        Eigen::VectorXd stim2;
        Eigen::VectorXd stim3;

    public:
        SAFEParams safe_params;

        Op_SAFE(GroptParams &_gparams, double _stim_thresh, double _weight_mod, bool _true_safe);
        virtual void init();

        virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void prox(Eigen::VectorXd &X);
        virtual void check(Eigen::VectorXd &X);

};

}  // close "namespace Gropt"

#endif
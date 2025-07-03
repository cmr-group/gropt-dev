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

    public:
        std::vector<double> tau1 = std::vector<double>(3, 0.0);
        std::vector<double> tau2 = std::vector<double>(3, 0.0);
        std::vector<double> tau3 = std::vector<double>(3, 0.0);
        std::vector<double> alpha1 = std::vector<double>(3, 0.0);
        std::vector<double> alpha2 = std::vector<double>(3, 0.0);
        std::vector<double> alpha3 = std::vector<double>(3, 0.0);
        std::vector<double> a1 = std::vector<double>(3, 0.0);
        std::vector<double> a2 = std::vector<double>(3, 0.0);
        std::vector<double> a3 = std::vector<double>(3, 0.0);
        std::vector<double> stim_limit = std::vector<double>(3, 0.0);
        std::vector<double> g_scale = std::vector<double>(3, 0.0);

        SAFEParams() = default;
};

class Op_SAFE : public Operator
{  
    protected:
        double stim_thresh;

        Eigen::VectorXd signs;
        Eigen::VectorXd stim1;
        Eigen::VectorXd stim2;
        Eigen::VectorXd stim3;

        SAFEParams safe_params;

    public:
        bool true_safe;

        Op_SAFE(GroptParams &_gparams, double _stim_thresh);
        virtual void init();

        void set_demo_params();

        virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void prox(Eigen::VectorXd &X);
        virtual void check(Eigen::VectorXd &X);

};

}  // close "namespace Gropt"

#endif
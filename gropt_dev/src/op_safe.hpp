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

class Op_SAFE : public Operator
{  
    protected:
        double stim_thresh;

        Eigen::VectorXd signs;
        Eigen::VectorXd stim1;
        Eigen::VectorXd stim2;
        Eigen::VectorXd stim3;

        Eigen::Vector3d tau1;
        Eigen::Vector3d tau2;
        Eigen::Vector3d tau3;
        Eigen::Vector3d alpha1;
        Eigen::Vector3d alpha2;
        Eigen::Vector3d alpha3;
        Eigen::Vector3d a1;
        Eigen::Vector3d a2;
        Eigen::Vector3d a3;
        Eigen::Vector3d stim_limit;
        Eigen::Vector3d g_scale;

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
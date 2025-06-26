#ifndef OP_BVALUE_H
#define OP_BVALUE_H

/**
 * Return identity matrix, to be used for regularization most likely.abort
 * i.e. simple duty cycle minimization can be accomplished with this.
 */

#include <iostream> 
#include <string>
#include <math.h>  
#include "Eigen/Dense"

#include "op_main.hpp"

namespace Gropt {

class Op_BValue : public Operator
{  
    protected:
        int start_ind = 0; 
        double GAMMA;
        double MAT_SCALE;

    public:
        Op_BValue(GroptParams &_gparams);
        virtual void init();

        virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void prox(Eigen::VectorXd &X);
        virtual void check(Eigen::VectorXd &X);
        double get_bvalue(Eigen::VectorXd &X);

};

}  // close "namespace Gropt"

#endif
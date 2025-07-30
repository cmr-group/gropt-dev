#ifndef OP_IDENTITY_H
#define OP_IDENTITY_H

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

class Op_Identity : public Operator
{  
    protected:

    public:
        Op_Identity(GroptParams &_gparams);
        Op_Identity(GroptParams &_gparams, double _weight_mod);
        
        virtual void init();

        virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void prox(Eigen::VectorXd &X);
        virtual void check(Eigen::VectorXd &X);

};

}  // close "namespace Gropt"

#endif
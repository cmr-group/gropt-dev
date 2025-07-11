#ifndef OP_SLEW_H
#define OP_SLEW_H

/**
 * Constriant on slew rate, i.e. |dG/dt| <= smax
 * Supports the 'rot_variant" option to constrain either the individual slews
 * or the slew magnitude.
 */

#include <iostream> 
#include <string>
#include <math.h>  
#include "Eigen/Dense"

#include "op_main.hpp"

namespace Gropt {

class Op_Slew : public Operator
{  
    protected:
        double smax;

    public:
        Op_Slew(GroptParams &_gparams, double _smax, bool _rot_variant, double _weight_mod);
        virtual void init();

        virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void prox(Eigen::VectorXd &X);
        virtual void check(Eigen::VectorXd &X);

};

}  // close "namespace Gropt"

#endif
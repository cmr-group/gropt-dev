#ifndef OP_TV_H
#define OP_TV_H

/**
 * Constriant on total variation, i.e. |dG|_1 <= TV
 * This is basically the slew rate constriant, but the prox function has been changed to 
 * a shrinkage operator.  They could probably be merged into one operator, or even considering
 * some way to make all constraints have shrinkage as a prox function option.
 */

#include <iostream> 
#include <string>
#include <math.h>  
#include "Eigen/Dense"

#include "op_main.hpp"

namespace Gropt {

class Op_TV : public Operator
{  
    protected:
        double weight_in = 1.0;;
        double tv_lam = 0.0;

    public:
        Op_TV(GroptParams &_gparams, double _tv_lam, double _weight_in);
        virtual void init();

        virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void prox(Eigen::VectorXd &X);
        virtual void check(Eigen::VectorXd &X);

};

}  // close "namespace Gropt"

#endif

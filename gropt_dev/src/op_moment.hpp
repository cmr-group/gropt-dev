#ifndef OP_MOMENT_H
#define OP_MOMENT_H

/**
 * Constraint on gradient moments, any order or moment
 * and different timing configurations and tolerances.
 */

#include <iostream> 
#include <string>
#include <vector>
#include "Eigen/Dense"

#include "op_main.hpp"

namespace Gropt {

class Op_Moment : public Operator
{  
    public:
        // These references may change due to gparams resize, so keep a record of input values
        int input_start = -1;
        int input_stop = -1;


        int moment_axis = 0;
        double moment_order;
        double moment_ref0 = 0.0;
        int moment_start;
        int moment_stop;
        double moment_target = 0;
        double moment_tol0 = 1e-6;

        Eigen::MatrixXd A;

        Op_Moment(GroptParams &_gparams, double _order);
        Op_Moment(GroptParams &_gparams, double _order, double _target);
        virtual void init();

        virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void prox(Eigen::VectorXd &X);
        virtual void check(Eigen::VectorXd &X);
};

}  // end namespace Gropt


#endif
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

    protected:
         // These references may change due to gparams resize, so keep a record of input values
        int start_idx0 = -1;
        int stop_idx0 = -1;
        int ref_idx0 = 0;

        int start_idx;
        int stop_idx;
        int ref_idx;    

        int moment_axis = 0;
        double moment_order;
        double moment_target = 0;
        double moment_tol0 = 1e-6;

        std::string units = "mT*ms/m";  

    public:
        Eigen::MatrixXd A;

        Op_Moment(GroptParams &_gparams, double _order, double _target, double _tol0, std::string _units, 
                  int _moment_axis, int _start_idx0, int _stop_idx0, int _ref_idx0, double _weight_mod);

        virtual void init();

        virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void prox(Eigen::VectorXd &X);
        virtual void check(Eigen::VectorXd &X);
};

}  // end namespace Gropt


#endif
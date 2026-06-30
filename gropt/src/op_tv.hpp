#ifndef OP_TV_H
#define OP_TV_H

#include <iostream>
#include <string>
#include <math.h>
#include "Eigen/Dense"

#include "op_main.hpp"

namespace Gropt {

class Op_TV : public Operator
{
    protected:
        double tv_lam = 0.0;
        int order = 1; // 1 = TV of the gradient (||slew||_1); 2 = TV of the slew (||jerk||_1)

    public:
        Op_TV(const ProblemData &_pdata, double _tv_lam, double _weight_mod, int _order = 1);
        virtual void init();

        virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void prox(Eigen::VectorXd &X);
        virtual void check(Eigen::VectorXd &X);

};

}  // close "namespace Gropt"

#endif

#ifndef OP_SLEW_H
#define OP_SLEW_H

/**
 * Slew-rate constraint |dG/dt| <= smax [T/m/s]: per axis (rot_variant) or on the cross-axis magnitude.
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
        Eigen::VectorXd smax_vec;

    public:
        Op_Slew(const ProblemData &_pdata, double _smax, bool _rot_variant, double _weight_mod);
        Op_Slew(const ProblemData &_pdata, const Eigen::VectorXd &_smax_vec, bool _rot_variant, double _weight_mod);
        virtual void init();

        virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void prox(Eigen::VectorXd &X);
        virtual void check(Eigen::VectorXd &X);

};

}  // close "namespace Gropt"

#endif

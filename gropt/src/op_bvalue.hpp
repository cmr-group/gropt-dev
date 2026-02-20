#ifndef OP_BVALUE_H
#define OP_BVALUE_H

#include <iostream>
#include <string>
#include <math.h>
#include "Eigen/Dense"

#include "op_main.hpp"

namespace Gropt {

enum BVALUE_MODE {
  BVALUE_MODE_SETVAL = 1,
  BVALUE_MODE_MINVAL = 2,
  BVALUE_MODE_MINVALMAX = 3,
};


class Op_BValue : public Operator
{
    protected:
        // These are the user inputted values
        int start_idx0 = -1;
        int stop_idx0 = -1;

        // These are the values after modifying problem parameters (i.e. changing N or Naxis)
        int start_idx;
        int stop_idx;

        // These are the values used in the calculations
        int i_start;
        int i_stop;

        double bval_target = 100;
        double bval_tol0 = 1;

        double GAMMA;
        double MAT_SCALE;

        BVALUE_MODE mode = BVALUE_MODE_SETVAL;
        double max_scale = 1.01;

    public:
        Op_BValue(const ProblemData &_pdata, double _bval_target, double _bval_tol0,
                  int _start_idx0, int _stop_idx0, double _weight_mod, BVALUE_MODE _mode, double _max_scale);
        virtual void init();

        virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void prox(Eigen::VectorXd &X);
        virtual void check(Eigen::VectorXd &X);
        double get_bvalue(Eigen::VectorXd &X);

};

}  // close "namespace Gropt"

#endif

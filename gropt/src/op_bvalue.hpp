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
        // Indices as passed by the caller (<= 0 = whole axis) and their working copies
        int start_idx0 = -1;
        int stop_idx0 = -1;

        int start_idx;
        int stop_idx;

        // Resolved [i_start, i_stop) range, set in init()
        int i_start;
        int i_stop;

        double bval_target = 100;
        double bval_tol0 = 1;

        double GAMMA;
        double MAT_SCALE;

        double PWL_A0, PWL_A1, PWL_A2; // pwl_quad: per-interval shifted-Legendre scales
        int ax_per_axis = 0;           // Ax length per axis: N, or 3(N-1) when pwl_quad

        BVALUE_MODE mode = BVALUE_MODE_SETVAL;
        double max_scale = 1.01;

    public:
        Op_BValue(const ProblemData &_pdata, double _bval_target, double _bval_tol0,
                  int _start_idx0, int _stop_idx0, double _weight_mod, BVALUE_MODE _mode, double _max_scale);
        virtual void init();

        // false (default): inclusive cumsum for q = integral g, then a rectangle sum for integral
        // q^2. Reads high by ~gamma^2 (dt^2/4) integral(g^2): 0.12% at 400 us, 0.0001% at 10 us.
        // true: exact for the piecewise-linear waveform the scanner plays, so b is raster
        // independent. q is piecewise quadratic, and on the shifted Legendre basis (orthogonal on
        // [0,dt], so no cross terms) with Q advanced by the trapezoid rule,
        //     c0 = Q_i + dt(2 g_i + g_i+1)/6   c1 = dt(g_i + g_i+1)/4   c2 = dt(g_i+1 - g_i)/12
        //     b  = gamma^2 dt sum_i [c0^2 + c1^2/3 + c2^2/5]
        // Still a perfect square, so b = ||Ax||^2 and prox/check/DCA are unchanged; Ax per axis
        // grows to 3(N-1) and spec_norm is estimated rather than guessed.
        // Off by default: it shifts every b-value target. Stability is a solver-settings question,
        // not a property of this flag -- on 24 hard PNS solves the bare default SolverCfg blew up
        // 1/24 with the rectangle rule and 3/24 with this one, a tuned recipe 0/24 with both.
        bool pwl_quad = false;

        virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void prox(Eigen::VectorXd &X);
        virtual void check(Eigen::VectorXd &X);
        double get_bvalue(Eigen::VectorXd &X);

};

}  // close "namespace Gropt"

#endif

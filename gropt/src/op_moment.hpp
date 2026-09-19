#ifndef OP_MOMENT_H
#define OP_MOMENT_H

/**
 * Constraint |M_k - target| <= tol on one gradient moment M_k = sum g(t) (t - t_ref)^k dt of any order k,
 * on one axis over an optional [start_idx, stop_idx) window. inv_vec applies the 180 sign flips.
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
        // Indices as passed by the caller (start/stop <= 0 = whole axis); init() uses the copies below
        int start_idx0 = -1;
        int stop_idx0 = -1;
        int ref_idx0 = 0;

        int start_idx;
        int stop_idx;
        int ref_idx;

        int moment_axis = 0;
        double moment_order;
        double moment_target = 0;
        double moment_tol0 = 1e-6;     // absolute tol, in this order's units
        double moment_tol0_m0 = 1e-6;  // tol in order-0 units, for the M0-anchored mode

        std::string units = "mT*ms/m";

    public:
        Eigen::MatrixXd A;

        // false (default): M0-anchored, tol is an order-0 tolerance scaled by ||A_k||/||A_0|| for order k.
        // true: tol is an absolute bound in this order's units.
        bool absolute_tol = false;

        Op_Moment(const ProblemData &_pdata, double _order, double _target, double _tol0, std::string _units,
                  int _moment_axis, int _start_idx0, int _stop_idx0, int _ref_idx0, double _weight_mod);

        virtual void init();

        virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void prox(Eigen::VectorXd &X);
        virtual void check(Eigen::VectorXd &X);

        // Appends the constant row A with the moment target (x0 is unused)
        virtual void append_eq_rows(std::vector<Eigen::VectorXd> &rows, std::vector<double> &targets,
                                    const Eigen::VectorXd &x0) const;
};

}  // end namespace Gropt


#endif

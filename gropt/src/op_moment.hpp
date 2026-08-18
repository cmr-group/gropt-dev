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
         // These references may change due to pdata resize, so keep a record of input values
        int start_idx0 = -1;
        int stop_idx0 = -1;
        int ref_idx0 = 0;

        int start_idx;
        int stop_idx;
        int ref_idx;

        int moment_axis = 0;
        double moment_order;
        double moment_target = 0;
        double moment_tol0 = 1e-6;     // absolute per-order tol (this order's physical units)
        double moment_tol0_m0 = 1e-6;  // M0-anchored reference tol (order-0 units; scaled by ||A_k||/||A_0|| in init)

        std::string units = "mT*ms/m";

    public:
        Eigen::MatrixXd A;

        // Tolerance mode. false (default) = M0-anchored: `tol` is the order-0 tolerance and higher orders
        // scale their bound by ||A_k||/||A_0||
        // true = the `tol` is an absolute bound in this order's physical units.
        bool absolute_tol = false;

        Op_Moment(const ProblemData &_pdata, double _order, double _target, double _tol0, std::string _units,
                  int _moment_axis, int _start_idx0, int _stop_idx0, int _ref_idx0, double _weight_mod);

        virtual void init();

        virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void prox(Eigen::VectorXd &X);
        virtual void check(Eigen::VectorXd &X);

        // Linear moment functional: a single constant row (A) with the moment target. eq_rows_vary
        // is false (the row does not depend on x0).
        virtual void append_eq_rows(std::vector<Eigen::VectorXd> &rows, std::vector<double> &targets,
                                    const Eigen::VectorXd &x0) const;
};

}  // end namespace Gropt


#endif

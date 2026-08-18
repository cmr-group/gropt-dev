#include "spdlog/spdlog.h"

#include "op_moment.hpp"

namespace Gropt {

Op_Moment::Op_Moment(const ProblemData &_pdata, double _order, double _target, double _tol0, std::string _units,
                     int _moment_axis, int _start_idx0, int _stop_idx0, int _ref_idx0, double _weight_mod)
    : Operator(_pdata) {
    name = "Moment";
    moment_order = _order;
    units = _units;

    double moment_scale = 1.0;
    if (units == "mT*ms/m") {
        moment_scale = 1.0;
    } else if (units == "T*s/m") {
        moment_scale = 1000.0 * pow(1000.0, moment_order + 1);
    } else if (units == "rad*s/m") {
        moment_scale = 1000.0 * pow(1000.0, moment_order + 1) / 4.257638544e7;
    } else if (units == "s/m") {
        moment_scale = 1000.0 * pow(1000.0, moment_order + 1) / 2.675153194e8;
    } else {
        spdlog::error("Unsupported units for moment constraint: {}", units);
        throw std::invalid_argument("Unsupported units for moment constraint");
    }

    // Order-0 units factor, for the M0-anchored tolerance (the default; see init). It uses the ORDER-0
    // units scale, not this order's, so the anchor doesn't double-scale the order (moment_scale already
    // carries the per-order factor for the absolute-tol mode and the target).
    double moment_scale0 = 1.0;
    if (units == "T*s/m") {
        moment_scale0 = 1000.0 * 1000.0;
    } else if (units == "rad*s/m") {
        moment_scale0 = 1000.0 * 1000.0 / 4.257638544e7;
    } else if (units == "s/m") {
        moment_scale0 = 1000.0 * 1000.0 / 2.675153194e8;
    } // mT*ms/m -> 1.0

    moment_target = _target * moment_scale;
    moment_tol0 = _tol0 * moment_scale;       // absolute per-order tol (this order's physical units)
    moment_tol0_m0 = _tol0 * moment_scale0;   // M0-anchored reference (order-0 units; scaled by ||A_k||/||A_0|| in init)
    moment_axis = _moment_axis;

    start_idx0 = _start_idx0;
    stop_idx0 = _stop_idx0;
    ref_idx0 = _ref_idx0;

    start_idx = start_idx0;
    stop_idx = stop_idx0;
    ref_idx = ref_idx0;

    weight_mod = _weight_mod;
}

void Op_Moment::init() {
    spdlog::trace("Op_Moment::init  N = {}", pdata->N);

    Ax_size = 1;

    A.setZero(1, pdata->Naxis * pdata->N);

    // If start and stop indices are not set, constraint covers the whole axis
    int i_start;
    if (start_idx <= 0) {
        i_start = moment_axis * pdata->N;
    } else {
        i_start = start_idx + moment_axis * pdata->N;
    }

    int i_stop;
    if (stop_idx <= 0) {
        i_stop = (moment_axis + 1) * pdata->N;
    } else {
        i_stop = stop_idx + moment_axis * pdata->N;
    }

    spec_norm2 = 0.0;
    for (int j = i_start; j < i_stop; j++) {
        double jj = j - moment_axis * pdata->N;
        double val = 1000.0 * 1000.0 * pdata->dt * pow((1000.0 * (pdata->dt * (jj - ref_idx))), moment_order);

        A(0, j) = val * pdata->inv_vec(j);
        spec_norm2 += val * val;
    }
    spec_norm = sqrt(spec_norm2);

    // Two tolerance modes (see absolute_tol):
    //  - M0-anchored (default): the passed tol is the ORDER-0 tolerance; higher orders scale their bound by
    //    ||A_k|| / ||A_0|| = spec_norm / spec_norm_0 (||A_0|| = the order-0 row norm over the same window,
    //    val_0 = 1e6*dt constant).
    //  - absolute_tol: the passed tol is taken directly in THIS order's physical units (moment_tol0 =
    //    tol * moment_scale(order)). Use it when you set a nonzero target for a higher order (e.g. a
    //    specified M2) and want the tolerance to mean an absolute physical bound, not a row-norm-scaled one.
    double base_val = 1000.0 * 1000.0 * pdata->dt; // val at order 0 (t^0 = 1)
    double spec_norm_0 = base_val * sqrt((double)(i_stop - i_start));

    target = moment_target;
    tol0 = absolute_tol ? moment_tol0 : moment_tol0_m0 * (spec_norm / spec_norm_0);
    tol = (1.0 - cushion) * tol0;

    if (do_init_weights) {
        obj_weight = 1.0;
        obj_weight *= weight_mod;
    }

    Operator::init();

    spdlog::trace("Initialized operator: {}", name);
    spdlog::trace("    moment_axis = {:d}  moment_order = {:.1f}", moment_axis, moment_order);
    spdlog::trace("    target = {:.1e}  tol0 = {:.1e}  tol = {:.1e}", target, tol0, tol);
    spdlog::trace("    i_start = {:d}  i_stop = {:d}", i_start, i_stop);
}

void Op_Moment::append_eq_rows(std::vector<Eigen::VectorXd> &rows, std::vector<double> &targets,
                               const Eigen::VectorXd &x0) const {
    if (!use_projection) return;
    rows.push_back(A.row(0).transpose()); // constant linear moment functional (independent of x0)
    targets.push_back(target);
}

void Op_Moment::forward(Eigen::VectorXd &X, Eigen::VectorXd &out) { out = A * X; }

void Op_Moment::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out) { out = A.transpose() * X; }

void Op_Moment::prox(Eigen::VectorXd &X) {
    spdlog::trace("Starting Op_Moment::prox");

    if (do_equil) {
        X.array() /= eq_rows.array();
    }
    X.array() *= spec_norm;

    for (int i = 0; i < X.size(); i++) {
        double lower_bound = (target - tol);
        double upper_bound = (target + tol);
        X(i) = X(i) < lower_bound ? lower_bound : X(i);
        X(i) = X(i) > upper_bound ? upper_bound : X(i);
    }

    if (do_equil) {
        X.array() *= eq_rows.array();
    }
    X.array() /= spec_norm;

    spdlog::trace("Finished Op_Moment::prox");
}

void Op_Moment::check(Eigen::VectorXd &X) {
    int is_feas = 1;

    if (do_equil) {
        X.array() /= eq_rows.array();
    }
    X.array() *= spec_norm;

    for (int i = 0; i < X.size(); i++) {
        if ((X(i) < target - tol0) || (X(i) > target + tol0)) {
            is_feas = 0;
        }
    }

    if (do_equil) {
        X.array() *= eq_rows.array();
    }
    X.array() /= spec_norm;

    hist_feas.push_back(is_feas);
}

} // namespace Gropt

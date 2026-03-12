#include "workspace_osqp.hpp"
#include "op_main.hpp"

namespace Gropt {

void WorkspaceOSQP::init(int Ax_size) {
    WorkspaceSolver::init(Ax_size);

    y00.setZero(Ax_size);
    z00.setZero(Ax_size);
}

void WorkspaceOSQP::reinit(int Ax_size) {
    WorkspaceSolver::reinit(Ax_size);

    if (y00.size() != Ax_size) {
        y00.setZero(Ax_size);
    }

    z00.setZero(Ax_size);
}

void WorkspaceOSQP::prep(Operator &op, const Eigen::VectorXd &X) {
    WorkspaceSolver::prep(op, X);
    z00 = z0;
}

void WorkspaceOSQP::init_reweight() {
    y00 = y0;
    z00 = z0;
}

void WorkspaceOSQP::reweight(int iiter, double scale_up, double scale_down, double eps) {
    double p = (y00 - y0).norm();
    double q = (z00 - z0).norm();

    if ((p < eps) && (q >= eps)) {
        weight_scale = 1 / scale_down;
    } else if ((p >= eps) && (q < eps)) {
        weight_scale = scale_up;
    } else if ((p >= eps) && (q >= eps)) {
        double new_weight = p / q;
        weight_scale = new_weight / weight;
    } else {
        weight_scale = 1.0;
    }

    // double max_exponent = std::log(1.0 + iiter);
    // max_exponent = pow(10.0, max_exponent);

    // if (weight > max_exponent) {
    //     weight = max_exponent;
    // } else if (weight < 1.0 / max_exponent) {
    //     weight = 1.0 / max_exponent;
    // }

    y00 = y0;
    z00 = z0;
}

} // namespace Gropt

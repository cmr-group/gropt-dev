#ifndef SDMM_WORKSPACE_H
#define SDMM_WORKSPACE_H

#include "Eigen/Dense"

namespace Gropt {

class Operator;  // Forward declaration

struct SDMMWorkspace {
    // Per-operator SDMM state
    double weight = 1.0;
    double gamma = 1.5;
    bool do_rw = true;
    bool do_gamma = true;
    bool do_weight = true;
    bool do_scalelim = true;

    Eigen::VectorXd y0, y1;
    Eigen::VectorXd z0, z1;
    Eigen::VectorXd s0, s1;
    Eigen::VectorXd yhat1, dyhat, dy, dhhat, dghat;
    Eigen::VectorXd yhat00, y00, s00, z00;

    void init(int Ax_size);
    void reinit(int Ax_size);
    void prep(Operator &op, const Eigen::VectorXd &X);
    void reweight(double rw_eps, double e_corr, double rw_scalelim);
};

}  // namespace Gropt

#endif

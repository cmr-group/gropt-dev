#include "workspace_sdmm.hpp"
#include "op_main.hpp"

namespace Gropt {

void WorkspaceSDMM::init(int Ax_size) {
    WorkspaceSolver::init(Ax_size);

    yhat1.setZero(Ax_size);
    dyhat.setZero(Ax_size);
    dy.setZero(Ax_size);
    dhhat.setZero(Ax_size);
    dghat.setZero(Ax_size);

    yhat00.setZero(Ax_size);
    y00.setZero(Ax_size);
    s00.setZero(Ax_size);
    z00.setZero(Ax_size);
}

void WorkspaceSDMM::reinit(int Ax_size) {
    WorkspaceSolver::reinit(Ax_size);

    if (yhat00.size() != Ax_size) {
        yhat00.setZero(Ax_size);
    }
    if (y00.size() != Ax_size) {
        y00.setZero(Ax_size);
    }

    yhat1.setZero(Ax_size);
    dyhat.setZero(Ax_size);
    dy.setZero(Ax_size);
    dhhat.setZero(Ax_size);
    dghat.setZero(Ax_size);

    s00.setZero(Ax_size);
    z00.setZero(Ax_size);
}

void WorkspaceSDMM::prep(Operator &op, const Eigen::VectorXd &X) {
    WorkspaceSolver::prep(op, X);
    z00 = z0;
}

void WorkspaceSDMM::reweight(double rw_eps, double e_corr, double rw_scalelim) {
    double rho0 = weight;

    yhat1.array() = y0.array() + rho0 * (s1.array() - z1.array());

    dyhat.array() = (yhat1.array() - yhat00.array());
    dy.array() = -(y1.array() - y00.array());
    dhhat.array() = s1.array() - s00.array();
    dghat.array() = -(z1.array() - z00.array());

    double norm_dhhat_dyhat = dhhat.norm() * dyhat.norm();
    double dot_dhhat_dhhat = dhhat.dot(dhhat);
    double dot_dhhat_dyhat = dhhat.dot(dyhat);

    double alpha_corr = 0.0;
    if ((norm_dhhat_dyhat > rw_eps) && (dot_dhhat_dhhat > rw_eps) && (dot_dhhat_dyhat > rw_eps)) {
        alpha_corr = dot_dhhat_dyhat / norm_dhhat_dyhat;
    }

    double norm_dghat_dy = dghat.norm() * dy.norm();
    double dot_dghat_dghat = dghat.dot(dghat);
    double dot_dghat_dy = dghat.dot(dy);

    double beta_corr = 0.0;
    if ((norm_dghat_dy > rw_eps) && (dot_dghat_dghat > rw_eps) && (dot_dghat_dy > rw_eps)) {
        beta_corr = dot_dghat_dy / norm_dghat_dy;
    }

    bool pass_alpha = false;
    bool pass_beta = false;

    double alpha = 0.0;
    if (alpha_corr > e_corr) {
        pass_alpha = true;
        double alpha_mg = dot_dhhat_dyhat / dot_dhhat_dhhat;
        double alpha_sd = dyhat.dot(dyhat) / dot_dhhat_dyhat;
        if (2.0 * alpha_mg > alpha_sd) {
            alpha = alpha_mg;
        } else {
            alpha = alpha_sd - 0.5 * alpha_mg;
        }
    }

    double beta = 0.0;
    if (beta_corr > e_corr) {
        pass_beta = true;
        double beta_mg = dot_dghat_dy / dot_dghat_dghat;
        double beta_sd = dy.dot(dy) / dot_dghat_dy;
        if (2.0 * beta_mg > beta_sd) {
            beta = beta_mg;
        } else {
            beta = beta_sd - 0.5 * beta_mg;
        }
    }

    double step_g1 = 0.0;
    double gamma1 = 0.0;
    if ((pass_alpha == true) && (pass_beta == true)) {
        step_g1 = sqrt(alpha * beta);
        gamma1 = 1.0 + 2.0 * sqrt(alpha * beta) / (alpha + beta);
    } else if ((pass_alpha == true) && (pass_beta == false)) {
        step_g1 = alpha;
        gamma1 = 1.9;
    } else if ((pass_alpha == false) && (pass_beta == true)) {
        step_g1 = beta;
        gamma1 = 1.1;
    } else {
        step_g1 = rho0;
        gamma1 = 1.5;
    }

    if (do_weight == true) {
        if ((do_scalelim == true) && (step_g1 > rw_scalelim * weight)) {
            weight *= rw_scalelim;
        } else if ((do_scalelim == true) && (rw_scalelim * step_g1 < weight)) {
            weight *= 1.0 / rw_scalelim;
        } else {
            weight = step_g1;
        }
    }

    if (do_gamma == true) {
        gamma = gamma1;
    }

    yhat00 = yhat1;
    y00 = y1;
    s00 = s1;
    z00 = z1;
}

} // namespace Gropt

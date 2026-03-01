#include "equilibrate.hpp"
#include "spdlog/spdlog.h"

namespace Gropt {

void equilibrate(GroptParams &gparams, int n_iter, int n_reps) {
    for (int i = 0; i < gparams.all_op.size(); i++) {
        Operator *op = gparams.all_op[i].get();

        op->eq_rows.setOnes(op->Ax_size);
        op->eq_cols.setOnes(gparams.N * gparams.Naxis);
        op->do_equil = true;
    }

    for (int iter = 0; iter < n_iter; iter++) {

        Eigen::VectorXd row_norms;
        Eigen::VectorXd col_norms;
        estimate_row_col_norms(gparams, n_reps, NormType::Inf, row_norms, col_norms);

        row_norms = row_norms.cwiseSqrt().cwiseInverse();
        col_norms = col_norms.cwiseSqrt().cwiseInverse();

        // TODO: This should be all axis
        col_norms[0] = col_norms[1];
        col_norms[col_norms.size() - 1] = col_norms[col_norms.size() - 2];

        int row_start = 0;
        for (int i = 0; i < gparams.all_op.size(); i++) {
            Operator *op = gparams.all_op[i].get();
            op->eq_rows.array() *= row_norms.segment(row_start, op->Ax_size).array();
            op->eq_cols.array() *= col_norms.array();
            row_start += op->Ax_size;
        }
    }
}

void get_eq_vecs(GroptParams &gparams, Eigen::VectorXd &row_norms, Eigen::VectorXd &col_norms) {

    int N_rows = 0;
    for (int i = 0; i < gparams.all_op.size(); i++) {
        N_rows += gparams.all_op[i]->Ax_size;
    }

    int N_cols = gparams.N * gparams.Naxis;

    row_norms.setZero(N_rows);
    col_norms.setZero(N_cols);

    int row_start = 0;
    for (int i = 0; i < gparams.all_op.size(); i++) {
        Operator *op = gparams.all_op[i].get();

        row_norms.segment(row_start, op->Ax_size) = op->eq_rows;
        if (i == 0) {
            col_norms = op->eq_cols;
        } else {
            // Check that col norms are the same across operators
            if (!col_norms.isApprox(op->eq_cols)) {
                spdlog::error("Column norms are not the same across operators!");
            }
        }
        row_start += op->Ax_size;
    }
}

void rescale_eq_vecs(GroptParams &gparams, double row_scale, double col_scale) {

    for (int i = 0; i < gparams.all_op.size(); i++) {
        Operator *op = gparams.all_op[i].get();
        op->eq_rows.array() *= row_scale;
        op->eq_cols.array() *= col_scale;
    }
}

void estimate_row_col_norms(GroptParams &gparams, int n_reps, NormType norm_type, Eigen::VectorXd &row_norms,
                            Eigen::VectorXd &col_norms) {
    spdlog::trace("Starting estimate_row_col_norms");

    // for (int i = 0; i < gparams.all_op.size(); i++) {
    //     Operator *op = gparams.all_op[i].get();
    //     spdlog::info("Operator {} has weight {} and spec_norm2 {}", op->name, op->obj_weight, op->spec_norm2);
    // }

    int N_rows = 0;
    for (int i = 0; i < gparams.all_op.size(); i++) {
        N_rows += gparams.all_op[i]->Ax_size;
    }

    int N_cols = gparams.N * gparams.Naxis;

    row_norms.setZero(N_rows);
    col_norms.setZero(N_cols);

    // Row norm reps
    for (int rep = 0; rep < n_reps; rep++) {
        Eigen::VectorXd z = Eigen::VectorXd::Random(N_cols).array().sign();
        int row_start = 0;
        for (int i = 0; i < gparams.all_op.size(); i++) {
            Operator *op = gparams.all_op[i].get();
            op->Ax_temp.setZero();
            op->forward_op(z, op->Ax_temp);
            if (norm_type == NormType::L2) {
                row_norms.segment(row_start, op->Ax_size) += op->Ax_temp.cwiseAbs2();
            } else {
                row_norms.segment(row_start, op->Ax_size) =
                    row_norms.segment(row_start, op->Ax_size).cwiseMax(op->Ax_temp.cwiseAbs());
            }
            row_start += op->Ax_size;
        }
    }

    // Col norm reps
    for (int rep = 0; rep < n_reps; rep++) {
        Eigen::VectorXd Atw;
        Atw.setZero(N_cols);
        for (int i = 0; i < gparams.all_op.size(); i++) {
            Operator *op = gparams.all_op[i].get();
            Eigen::VectorXd w = Eigen::VectorXd::Random(op->Ax_size).array().sign();
            op->x_temp.setZero();
            op->transpose_op(w, op->x_temp, false);
            Atw.array() += op->x_temp.array();
        }

        if (norm_type == NormType::L2) {
            col_norms += Atw.cwiseAbs2();
        } else {
            col_norms = col_norms.cwiseMax(Atw.cwiseAbs());
        }
    }

    if (norm_type == NormType::L2) {
        row_norms = (row_norms / n_reps).cwiseSqrt();
        col_norms = (col_norms / n_reps).cwiseSqrt();
    }

    spdlog::trace("Finished estimate_row_col_norms");
}

double estimate_spec_norm(GroptParams &gparams, int n_iters) {

    int Ntot = gparams.N * gparams.Naxis;

    // Make a VectorXd of size Ntot using a standard normal distribution, and normalize it to have norm 1
    Eigen::VectorXd v = Eigen::VectorXd::Random(Ntot);
    v.normalize();

    double lambda_max;
    for (int iter = 0; iter < n_iters; iter++) {
        Eigen::VectorXd v_next = Eigen::VectorXd::Zero(Ntot);
        for (int i = 0; i < gparams.all_op.size(); i++) {
            Operator *op = gparams.all_op[i].get();

            op->Ax_temp.setZero();
            op->x_temp.setZero();
            op->forward_op(v, op->Ax_temp);
            op->transpose_op(op->Ax_temp, op->x_temp);
            v_next += op->x_temp;
        }
        lambda_max = v_next.norm();
        v = v_next / lambda_max;
        spdlog::debug("estimate_spec_norm iteration {}: lambda_max = {}", iter, lambda_max);
    }

    return sqrt(lambda_max);
}
} // namespace Gropt

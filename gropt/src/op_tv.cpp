#include "spdlog/spdlog.h"

#include "op_tv.hpp"

namespace Gropt {

Op_TV::Op_TV(const ProblemData &_pdata, double _tv_lam, double _weight_mod, int _order)
    : Operator(_pdata)
{
    name = "TotalVariation";
    weight_mod = _weight_mod;
    tv_lam = _tv_lam;
    order = (_order == 2) ? 2 : 1;
}

void Op_TV::init()
{
    spdlog::trace("Op_TV::init  N = {}  order = {}", pdata->N, order);

    target = 0;
    tol0 = tv_lam;
    tol = (1.0-cushion) * tol0;

    if (order == 2) {
        // Second difference (jerk), stencil [1,-2,1]/dt^2: ||stencil|| = 4 -> operator norm 4/dt^2.
        spec_norm2 = 16.0 / (pdata->dt*pdata->dt*pdata->dt*pdata->dt);
        Ax_size = pdata->Naxis * (pdata->N-2);
    } else {
        // First difference (slew), stencil [-1,1]/dt: ||stencil|| = 2 -> operator norm 2/dt.
        spec_norm2 = 4.0/pdata->dt/pdata->dt;
        Ax_size = pdata->Naxis * (pdata->N-1);
    }
    spec_norm = sqrt(spec_norm2);

    if (do_init_weights) {
        obj_weight = 1.0;
        obj_weight *= weight_mod;
    }

    Operator::init();
}

void Op_TV::forward(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    if (order == 2) {
        // jerk = second difference of the gradient
        for (int i_ax = 0; i_ax < Naxis; i_ax++) {
            for (int i = 0; i < (N-2); i++) {
                out(i_ax*(N-2)+i) =
                    (X(i_ax*N+i+2) - 2.0*X(i_ax*N+i+1) + X(i_ax*N+i)) / (pdata->dt*pdata->dt);
            }
        }
    } else {
        // slew = first difference of the gradient
        for (int i_ax = 0; i_ax < Naxis; i_ax++) {
            for (int i = 0; i < (N-1); i++) {
                out(i_ax*(N-1)+i) = (X(i_ax*N+i+1) - X(i_ax*N+i))/pdata->dt;
            }
        }
    }
}

void Op_TV::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    if (order == 2) {
        // Adjoint of the second-difference operator. Row i (i in 0..N-3) has [1,-2,1]/dt^2 at
        // columns i,i+1,i+2; so out_j gets X_j - 2 X_{j-1} + X_{j-2} over the valid row range.
        double idt2 = 1.0 / (pdata->dt*pdata->dt);
        int M = N - 2; // rows per axis
        for (int i_ax = 0; i_ax < Naxis; i_ax++) {
            int xb = i_ax*N;
            int yb = i_ax*M;
            for (int j = 0; j < N; j++) {
                double term = 0.0;
                if (j <= M-1)                  term += X(yb + j);
                if ((j-1 >= 0) && (j-1 <= M-1)) term += -2.0 * X(yb + j-1);
                if ((j-2 >= 0) && (j-2 <= M-1)) term += X(yb + j-2);
                out(xb + j) = term * idt2;
            }
        }
    } else {
        for (int i_ax = 0; i_ax < Naxis; i_ax++) {
            out(i_ax*N+0) = -X(i_ax*(N-1)+0) / pdata->dt;
            for (int i = 1; i < (N-1); i++) {
                out(i_ax*N+i) = (X(i_ax*(N-1)+i-1) - X(i_ax*(N-1)+i)) / pdata->dt;
            }
            out(i_ax*N+N-1) = X(i_ax*(N-1)+N-2) / pdata->dt;
        }
    }
}

void Op_TV::prox(Eigen::VectorXd &X)
{
    spdlog::trace("Starting Op_TV::prox");

    // Prox of the L1 PENALTY  f(u) = tv_lam * ||u||_1  on the physical slew (order 1) / jerk
    // (order 2)
    double thresh = tv_lam * spec_norm2 / rho;

    if (do_equil) {
        X.array() /= eq_rows.array();
    }
    X.array() *= spec_norm;

    for (int i = 0; i < X.size(); i++) {
        if (fabs(X(i)) > thresh) {
            X(i) = X(i) > 0 ? X(i) - thresh : X(i) + thresh;
        } else {
            X(i) = 0.0;
        }
    }

    if (do_equil) {
        X.array() *= eq_rows.array();
    }
    X.array() /= spec_norm;

    spdlog::trace("Finished Op_TV::prox");
}


void Op_TV::check(Eigen::VectorXd &X)
{
    int is_feas = 1;

    // As of right now we will assume that the TV operator is always feasible
    // This isn't necessarily the case, we should have an option to check feasibility (TODO)

    hist_feas.push_back(is_feas);
}

}  // close "namespace Gropt"

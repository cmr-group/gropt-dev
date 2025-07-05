#include "spdlog/spdlog.h"

#include "op_safe.hpp"

namespace Gropt {

Op_SAFE::Op_SAFE(GroptParams &_gparams, double _stim_thresh)
    : Operator(_gparams)
{
    name = "SAFE"; 
    stim_thresh = _stim_thresh;

    // Proper SAFE model is not stable in Krylov optimization, Setting to false removes the abs(), which
    // makes it converge nicely, but does not give the proper SAFE, though in every example I have tried
    // it still matches *PEAK* SAFE (it is the lower safe values that might not match up)
    true_safe = false; 
}

void Op_SAFE::init()
{
    spdlog::trace("Op_SAFE::init  N = {}", gparams->N);

    safe_params.calc_alphas(gparams->dt);

    target = 0;
    tol0 = stim_thresh;
    tol = (1.0-cushion) * tol0;

    spec_norm2 = 4.0/gparams->dt/gparams->dt;
    spec_norm = sqrt(spec_norm2);

    Ax_size = 3 * gparams->Naxis * gparams->N;

    if (!save_weights) {
        weight = 1e4 / spec_norm2;
    }

    signs.setZero(gparams->Naxis*gparams->N);
    stim1.setZero(gparams->Naxis*gparams->N);
    stim2.setZero(gparams->Naxis*gparams->N);
    stim3.setZero(gparams->Naxis*gparams->N);

    Operator::init();

    spdlog::trace("Op_SAFE::init  Done!");
}


void SAFEParams::calc_alphas(double dt) 
{
    for (int i = 0; i < 3; i++) {
        alpha1[i] = dt/(tau1[i] + dt);
        alpha2[i] = dt/(tau2[i] + dt);
        alpha3[i] = dt/(tau3[i] + dt);
    }
}


void Op_SAFE::forward(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    for (int j = 0; j < Naxis; j++) {
        out(j*N) = X(j*N)/dt;
        for (int i = 1; i < N; i++) {
            out(j*N+i) = (X(j*N+i) - X(j*N+i-1))/dt;
        }
    }

    for (int i = 0; i < signs.size(); i++) {
        if (out(i) < 0) {signs(i) = -1.0;}
        else {signs(i) = 1.0;}
    }

    stim1.setZero();
    for (int j = 0; j < Naxis; j++) {
        stim1(j*N) = safe_params.alpha1[j] * out(j*N);
        for (int i = 1; i < N; i++) {
            stim1(j*N+i) = safe_params.alpha1[j] * out(j*N+i) + (1.0-safe_params.alpha1[j]) * stim1(j*N+i-1);
        }
    }
    if (true_safe) {
        for (int i = 0; i < stim1.size(); i++) {
            stim1(i) = abs(stim1(i));
        }
    }

    stim2.setZero();
    for (int j = 0; j < Naxis; j++) {
        stim2(j*N) = safe_params.alpha2[j] * abs(out(j*N));
        for (int i = 1; i < N; i++) {
            stim2(j*N+i) = safe_params.alpha2[j] * abs(out(j*N+i)) + (1.0-safe_params.alpha2[j]) * stim2(j*N+i-1);
        }
    }
    if (!true_safe) {
        for (int i = 0; i < stim2.size(); i++) {
            stim2(i) = signs(i) * stim2(i);
        }
    }

    stim3.setZero();
    for (int j = 0; j < Naxis; j++) {
        stim3(j*N) = safe_params.alpha3[j] * out(j*N);
        for (int i = 1; i < N; i++) {
            stim3(j*N+i) = safe_params.alpha3[j] * out(j*N+i) + (1.0-safe_params.alpha3[j]) * stim3(j*N+i-1);
        }
    }
    if (true_safe) {
        for (int i = 0; i < stim3.size(); i++) {
            stim3(i) = abs(stim3(i)); 
        }
    }

    // for (int j = 0; j < Naxis; j++) {
    //     for (int i = 0; i < N; i++) {
    //         out(j*N+i) = (a1(j)*stim1(j*N+i) + a2(j)*stim2(j*N+i) + a3(j)*stim3(j*N+i)) / stim_limit(j) * g_scale(j);
    //     }
    // }

    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            out(j*3*N+i) = safe_params.a1[j] * stim1(j*N+i) / safe_params.stim_limit[j] * safe_params.g_scale[j];
            out(j*3*N+i+N) = safe_params.a2[j] * stim2(j*N+i) / safe_params.stim_limit[j] * safe_params.g_scale[j];
            out(j*3*N+i+2*N) = safe_params.a3[j] * stim3(j*N+i) / safe_params.stim_limit[j] * safe_params.g_scale[j];
        }
    }

}

void Op_SAFE::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out)
{
    x_temp.setZero();
    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            x_temp(j*N+i) = X(j*3*N+i) + X(j*3*N+i+N) + X(j*3*N+i+2*N);
        }
    }

    
    for (int i = 0; i < signs.size(); i++) {
        if (x_temp(i) < 0) {signs(i) = -1.0;}
        else {signs(i) = 1.0;}
    }

    // for (int j = 0; j < Naxis; j++) {
    //     for (int i = 0; i < N; i++) {
    //         if (X(j*3*N+i+N) < 0) {signs(j*N+i) = -1.0;}
    //         else {signs(j*N+i) = 1.0;}
    //     }
    // }


    stim1.setZero();
    for (int j = 0; j < Naxis; j++) {
        stim1(j*N+N-1) = safe_params.alpha1[j] * x_temp(j*N+N-1);
        for (int i = N-2; i >= 0; i--) {
            stim1(j*N+i) = safe_params.alpha1[j] * x_temp(j*N+i) + (1-safe_params.alpha1[j]) * stim1(j*N+i+1);
        }
    }
    if (true_safe) {
        for (int i = 0; i < x_temp.size(); i++) {
            stim1(i) = abs(stim1(i));
        }
    }


    stim2.setZero();
    for (int j = 0; j < Naxis; j++) {
        stim2(j*N+N-1) = safe_params.alpha2[j] * abs(x_temp(j*N+N-1));
        for (int i = N-2; i >= 0; i--) {
            stim2(j*N+i) = safe_params.alpha2[j] * abs(x_temp(j*N+i)) + (1-safe_params.alpha2[j]) * stim2(j*N+i+1);
        }
    }
    if (!true_safe) {
        for (int i = 0; i < x_temp.size(); i++) {
            stim2(i) = signs(i) * stim2(i);
        }
    }



    stim3.setZero();
    for (int j = 0; j < Naxis; j++) {
        stim3(j*N+N-1) = safe_params.alpha3[j] * x_temp(j*N+N-1);
        for (int i = N-2; i >= 0; i--) {
            stim3(j*N+i) = safe_params.alpha3[j] * x_temp(j*N+i) + (1-safe_params.alpha3[j]) * stim3(j*N+i+1);
        }
    }
    if (true_safe) {
        for (int i = 0; i < x_temp.size(); i++) {
            stim3(i) = abs(stim3(i));
        }
    }


    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            // out(j*N+i) = signs(j*N+i) * (a1(j)*stim1(j*N+i) + a2(j)*stim2(j*N+i) + a3(j)*stim3(j*N+i)) / stim_limit(j) * g_scale(j);
            out(j*N+i) = (safe_params.a1[j]*stim1(j*N+i) + 
                          safe_params.a2[j]*stim2(j*N+i) + 
                          safe_params.a3[j]*stim3(j*N+i)) / safe_params.stim_limit[j] * safe_params.g_scale[j];
        }
    }


    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N-1; i++) {
            out(j*N+i) = (out(j*N+i) - out(j*N+i+1))/dt;
        }
        out(j*N+N-1) = out(j*N+N-1)/dt;
    }
}

void Op_SAFE::prox(Eigen::VectorXd &X)
{
    spdlog::trace("Starting Op_SAFE::prox");

    // This just sums the three terms, it is not the cross axis sum.
    x_temp.setZero();
    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            // x_temp(j*N+i) = X(j*3*N+i) + X(j*3*N+i+N) + X(j*3*N+i+2*N);
            x_temp(j*N+i) = abs(X(j*3*N+i)) + abs(X(j*3*N+i+N)) + abs(X(j*3*N+i+2*N));
        }
    }

    double upper_bound = (target+tol);

    if (rot_variant) {
        for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            double val = abs(x_temp(j*N+i));

            if (val > upper_bound) {
                X(j*3*N+i) *= (upper_bound/val);
                X(j*3*N+i+N) *= (upper_bound/val);
                X(j*3*N+i+2*N) *= (upper_bound/val);
            }
        }
        }   
    } else {
        for (int i = 0; i < N; i++) {
            double val = 0.0;
            for (int j = 0; j < Naxis; j++) {
                val += X(j*N+i)*X(j*N+i);
            }
            val = sqrt(val);

            if (val > upper_bound) {
                for (int j = 0; j < Naxis; j++) {
                    X(j*N+i) *= (upper_bound/val);
                }
            }
        }
    }

    spdlog::trace("Finished Op_SAFE::prox");
}


void Op_SAFE::check(Eigen::VectorXd &X)
{
    int is_feas = 1;

    x_temp.setZero();
    for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            // x_temp(j*N+i) = X(j*3*N+i) + X(j*3*N+i+N) + X(j*3*N+i+2*N);
            x_temp(j*N+i) = abs(X(j*3*N+i)) + abs(X(j*3*N+i+N)) + abs(X(j*3*N+i+2*N));
        }
    }

    double upper_bound = (target+tol0);

    if (rot_variant) {
        for (int j = 0; j < Naxis; j++) {
        for (int i = 0; i < N; i++) {
            double val = abs(x_temp(j*N+i));

            if (val > upper_bound) {
                is_feas = 0;
            }
        }
        }   
    } else {
        for (int i = 0; i < N; i++) {
            double val = 0.0;
            for (int j = 0; j < Naxis; j++) {
                val += X(j*N+i)*X(j*N+i);
            }
            val = sqrt(val);

            if (val > upper_bound) {
                is_feas = 0;
            }
        }
    }

    hist_feas.push_back(is_feas);
}


void SAFEParams::set_demo_params() 
{
    tau1[0] = 0.135/1000.0; 
    tau2[0] = 12.0/1000.0;
    tau3[0] = 0.5175/1000.0;
    a1[0] = 0.3130;
    a2[0] = 0.1995;
    a3[0] = 0.4875;
    stim_limit[0] = 27.894;
    g_scale[0] = 0.325;

    tau1[1] = 1.5/1000.0;
    tau2[1] = 2.5/1000.0;
    tau3[1] = 0.15/1000.0;
    a1[1] = 0.55;
    a2[1] = 0.15;
    a3[1] = 0.3;
    stim_limit[1] = 15;
    g_scale[1] = 0.31;

    tau1[2] = 2/1000.0;
    tau2[2] = 0.12/1000.0;
    tau3[2] = 1/1000.0;
    a1[2] = 0.42;
    a2[2] = 0.4;
    a3[2] = 0.18;
    stim_limit[2] = 25;
    g_scale[2] = 0.25;
}

void SAFEParams::set_params(double *_tau1, double *_tau2, double *_tau3, 
                            double *_a1, double *_a2, double *_a3,
                            double *_stim_limit, double *_g_scale) 
{
    tau1[0] = _tau1[0]; 
    tau2[0] = _tau2[0];
    tau3[0] = _tau3[0];
    a1[0] = _a1[0];
    a2[0] = _a2[0];
    a3[0] = _a3[0];
    stim_limit[0] = _stim_limit[0];
    g_scale[0] = _g_scale[0];

    tau1[1] = _tau1[1];
    tau2[1] = _tau2[1];
    tau3[1] = _tau3[1];
    a1[1] = _a1[1];
    a2[1] = _a2[1];
    a3[1] = _a3[1];
    stim_limit[1] = _stim_limit[1];
    g_scale[1] = _g_scale[1];

    tau1[2] = _tau1[2];
    tau2[2] = _tau2[2];
    tau3[2] = _tau3[2];
    a1[2] = _a1[2];
    a2[2] = _a2[2];
    a3[2] = _a3[2];
    stim_limit[2] = _stim_limit[2];
    g_scale[2] = _g_scale[2];
}

void SAFEParams::swap_first_axes(int new_first_axis) 
{
    if (new_first_axis < 0 || new_first_axis >= 3) {
        spdlog::error("Invalid axis index: {}", new_first_axis);
        return;
    } else if (new_first_axis > 0) {
        std::swap(tau1[0], tau1[new_first_axis]);
        std::swap(tau2[0], tau2[new_first_axis]);
        std::swap(tau3[0], tau3[new_first_axis]);
        std::swap(a1[0], a1[new_first_axis]);
        std::swap(a2[0], a2[new_first_axis]);
        std::swap(a3[0], a3[new_first_axis]);
        std::swap(stim_limit[0], stim_limit[new_first_axis]);
        std::swap(g_scale[0], g_scale[new_first_axis]);
    }
}

}  // close "namespace Gropt"
#include "spdlog/spdlog.h"

#include "op_main.hpp"

namespace Gropt {

    Operator::Operator(GroptParams &_gparams)
    {
        name = "OperatorMain";
        gparams = &_gparams;
    } 

    Operator::~Operator() {}

    // TODO: This also needs to set things like weight and gamma to their initial values
    void Operator::init()
    {
        spdlog::trace("In Operator::init() from {}", name);

        N = gparams->N;
        Naxis = gparams->Naxis;
        dt = gparams->dt;
        Ntot = N * Naxis;

        x_temp.setZero(Ntot);
        x_temp_obj.setZero(Ntot);
        Ax_temp.setZero(Ax_size);
    }

    void Operator::init_parsdmm()
    {
        y0.setZero(Ax_size);
        y1.setZero(Ax_size);
        z0.setZero(Ax_size);
        z1.setZero(Ax_size);

        s0.setZero(Ax_size);
        s1.setZero(Ax_size);

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

    void Operator::reinit_parsdmm()
    {
        // y0.setZero(Ax_size);
        // y1.setZero(Ax_size);

        s0.setZero(Ax_size);
        s1.setZero(Ax_size);

        yhat1.setZero(Ax_size);
        dyhat.setZero(Ax_size);
        dy.setZero(Ax_size);
        dhhat.setZero(Ax_size);
        dghat.setZero(Ax_size);

        // yhat00.setZero(Ax_size);
        s00.setZero(Ax_size);
        z00.setZero(Ax_size);
    }

    void Operator::prep_parsdmm(Eigen::VectorXd &X)
    {
        forward_op(X, z0);
        z1 = z0;
        z00 = z0;
    }


    void Operator::forward(Eigen::VectorXd &X, Eigen::VectorXd &out)
    {
        spdlog::warn("Operator::forward is not implemented for the base Operator class.");
    }

    void Operator::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out)
    {
        spdlog::warn("Operator::transpose is not implemented for the base Operator class.");
    }

    void Operator::forward_op(Eigen::VectorXd &X, Eigen::VectorXd &out)
    {
        forward(X, out);
    }

    void Operator::transpose_op(Eigen::VectorXd &X, Eigen::VectorXd &out)
    {
        transpose(X, out);
        out.array() *= gparams->fixer.array();
        // out.array() /= spec_norm2;
    }

    void Operator::add_Atb(Eigen::VectorXd &b)
    {
        spdlog::trace("Operator::add_Atb  start  name = {}", name);
        
        Ax_temp.setZero();
        x_temp.setZero();

        if (gparams->solver_method == SolverMethod::GROPT_SDMM) {
            Ax_temp = weight*z0 - y0;
        }
    
        transpose_op(Ax_temp, x_temp);
        b += x_temp;

        spdlog::trace("Operator::add_Atb  end    name = {}", name);
    }

    void Operator::add_AtAx(Eigen::VectorXd &X, Eigen::VectorXd &out)
    {
        spdlog::trace("Operator::add_AtAx  start  name = {}", name);
        
        Ax_temp.setZero();
        x_temp.setZero();
        forward_op(X, Ax_temp);
        transpose_op(Ax_temp, x_temp);

        out.array() += weight*x_temp.array();

        spdlog::trace("Operator::add_AtAx  end    name = {}", name);
    }

    void Operator::add_obj(Eigen::VectorXd &X, Eigen::VectorXd &out)
    {
        Ax_temp.setZero();
        x_temp.setZero();

        forward_op(X, Ax_temp);
        transpose_op(Ax_temp, x_temp);

        out.array() += obj_weight*Ax_temp.array();

        spdlog::trace("Operator::add_obj   name = {}  obj_weight = {:.1e}", name, obj_weight);
    }

    void Operator::check(Eigen::VectorXd &X)
    {
        feas_check = (X.array() - target).abs().maxCoeff();

        if (feas_check <= tol0) {
            hist_feas.push_back(1);
        } else {
            hist_feas.push_back(0);
        }
    }

    void Operator::prox(Eigen::VectorXd &X)
    {
        spdlog::warn("Operator::prox is not implemented for the base Operator class.");
    }

    void Operator::get_feas(Eigen::VectorXd &s)
    {
        feas_temp = s;
        prox(feas_temp);
        feas_temp = s - feas_temp;

        r_feas = feas_temp.cwiseAbs().maxCoeff()/(s.cwiseAbs().maxCoeff() + 1.0e-32);

        hist_r_feas.push_back(r_feas);
    }

    void Operator::reweight_parsdmm()
    {
        double rho0 = weight;
        
        // Everthing with a "y" needs to be double checked for negation
        yhat1.array() = y0.array() + rho0*(s1.array() - z1.array());
        
        dyhat.array() = (yhat1.array() - yhat00.array());
        dy.array() = -(y1.array() - y00.array());
        dhhat.array() = s1.array() - s00.array();
        dghat.array() = -(z1.array() - z00.array());


        double norm_dhhat_dyhat = dhhat.norm()*dyhat.norm();
        double dot_dhhat_dhhat = dhhat.dot(dhhat);
        double dot_dhhat_dyhat = dhhat.dot(dyhat);

        
        double alpha_corr = 0.0;
        if ((norm_dhhat_dyhat > rw_eps) 
            && (dot_dhhat_dhhat > rw_eps) 
            && (dot_dhhat_dyhat > rw_eps)) {
                alpha_corr = dot_dhhat_dyhat/norm_dhhat_dyhat;
            }

        double norm_dghat_dy = dghat.norm()*dy.norm();
        double dot_dghat_dghat = dghat.dot(dghat);
        double dot_dghat_dy = dghat.dot(dy);

        double beta_corr = 0.0;
        if ((norm_dghat_dy > rw_eps) 
            && (dot_dghat_dghat > rw_eps) 
            && (dot_dghat_dy > rw_eps)) {
                beta_corr = dot_dghat_dy/norm_dghat_dy;
            }

        bool pass_alpha = false;
        bool pass_beta = false;

        double alpha = 0.0;
        if (alpha_corr > e_corr) {
            pass_alpha = true;
            double alpha_mg = dot_dhhat_dyhat/dot_dhhat_dhhat;
            double alpha_sd = dyhat.dot(dyhat)/dot_dhhat_dyhat;
            if (2.0*alpha_mg > alpha_sd) {
                alpha = alpha_mg;
            } else {
                alpha = alpha_sd - 0.5*alpha_mg;
            }
        }

        double beta = 0.0;
        if (beta_corr > e_corr) {
            pass_beta = true;
            double beta_mg = dot_dghat_dy/dot_dghat_dghat;
            double beta_sd = dy.dot(dy)/dot_dghat_dy;
            if (2.0*beta_mg > beta_sd) {
                beta = beta_mg;
            } else {
                beta = beta_sd - 0.5*beta_mg;
            }
        }

        double step_g1 = 0.0;
        double gamma1 = 0.0;
        if ((pass_alpha == true) && (pass_beta == true)) {
            step_g1 = sqrt(alpha*beta);
            gamma1 = 1.0 + 2.0*sqrt(alpha*beta)/(alpha+beta);
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
            if ((do_scalelim == true) && (step_g1 > rw_scalelim*weight)) {
                weight *= rw_scalelim;
            } else if ((do_scalelim == true) && (rw_scalelim*step_g1 < weight)) {
                weight *= 1.0/rw_scalelim;
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


}  // close "namespace Gropt"
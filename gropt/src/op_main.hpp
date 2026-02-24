#ifndef OP_MAIN_H
#define OP_MAIN_H

/**
 * This is our main parent class for every constraint and regularization term in GrOpt
 * The main functions that we use here are the reweighting schemes, where most operators will
 * then just need to implement the forward, transpose, and proximal mapping operators.
 */

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include "Eigen/Dense"

#include "problem_data.hpp"

namespace Gropt {

struct WorkspaceSolver;  // Forward declaration

class Operator  // This is the main parent class for every operator in GrOpt
{
    public:
        std::string name;

        const ProblemData *pdata;

        int N;
        int Naxis;
        int Ntot;
        double dt;
        int Ax_size;

        bool rot_variant = true;
        bool do_init_weights = true;

        double target;
        double tol0;
        double tol;
        double cushion = 1e-2;  // This is the cushion factor, which is used to reduce the tolerance by a factor of 1-cushion

        double spec_norm2;
        double spec_norm;

        double obj_weight = 1.0;
        double weight_mod = 1.0;

        Eigen::VectorXd x_temp;
        Eigen::VectorXd x_temp_obj;
        Eigen::VectorXd Ax_temp;

        double feas_check;
        double r_feas;
        Eigen::VectorXd feas_temp;

        std::vector<int> hist_feas;
        std::vector<double> hist_r_feas;

        // ----------------------------------------
        Operator() = default;
        Operator(const ProblemData &_pdata);
        virtual ~Operator();

        virtual void init();

        virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void forward_op(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose_op(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void check(Eigen::VectorXd &X);
        virtual void prox(Eigen::VectorXd &X);
        virtual void get_feas(Eigen::VectorXd &s);
        virtual void add_Atb(Eigen::VectorXd &b, const WorkspaceSolver &ws);
        void add_AtAx(Eigen::VectorXd &x, Eigen::VectorXd &out, const WorkspaceSolver &ws);
        virtual void add_obj(Eigen::VectorXd &x, Eigen::VectorXd &out);
};

}  // close "namespace Gropt"

#endif

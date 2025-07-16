#include "spdlog/spdlog.h"
#include <matplot/matplot.h> 
#include <highfive/H5Easy.hpp>

#include <vector>
#include <algorithm>

#include "gropt_params.hpp"
#include "solver.hpp"
#include "solver_groptsdmm.hpp"
#include "op_gradient.hpp"
#include "op_slew.hpp"
#include "op_moment.hpp"
#include "op_identity.hpp"
#include "op_bvalue.hpp"


using namespace Gropt;

void save_debug(GroptParams &params, Solver &solver) {
    spdlog::info("Saving debug information...");

    H5Easy::File file("debug_output.h5", H5Easy::File::Overwrite);
    H5Easy::dump(file, "/hist_x", solver.hist_X);

    H5Easy::dump(file, "/hist_cg_iter", solver.hist_cg_iter);

    // for (int i = 0; i < params.all_op.size(); i++) {
    //     Operator *op = params.all_op[i];
    //     H5Easy::dump(file, "/op/" + op->name + "/hist_Ax", op->hist_Ax);
    //     H5Easy::dump(file, "/op/" + op->name + "/hist_y", op->hist_y);
    //     H5Easy::dump(file, "/op/" + op->name + "/hist_z", op->hist_z);
    //     H5Easy::dump(file, "/op/" + op->name + "/hist_Aty", op->hist_Aty);
    // }

    // for (int i = 0; i < params.all_obj.size(); i++) {
    //     Operator *op = params.all_obj[i];
    //     H5Easy::dump(file, "/obj/" + op->name + "/hist_obj", op->hist_obj);
    // }

    spdlog::info("Debug information saved.");
}

/*
void simple_test() {
    GroptParams params;
    params.N = 103;
    params.Naxis = 1;
    params.dt = 10e-6;
    params.vec_init_simple();

    double gmax = 30e-3;
    params.all_op.push_back(new Op_Gradient(params, gmax));

    double smax = 200;
    params.all_op.push_back(new Op_Slew(params, smax));

    params.all_op.push_back(new Op_Moment(params, 0, 2.0));
    params.all_op.push_back(new Op_Moment(params, 1, 0.0));
    params.all_op.push_back(new Op_Moment(params, 2, 0.0));
    params.all_op.push_back(new Op_Moment(params, 3, 0.0));  // N = 103 for feasibility

    // params.all_obj.push_back(new Op_Identity(params));  

    // SolverPARSDMM solver(params);
    SolverGroptSDMM solver(params);
    solver.N_iter = 300;

    params.init();
    solver.solve();

    std::vector<double> vec(params.final_X.data(), params.final_X.data() + params.final_X.size());
    matplot::plot(vec);
    matplot::show();

    spdlog::info("Finished main");
}

void maximize_I() 
{
    GroptParams params;
    params.N = 103;
    params.Naxis = 1;
    params.dt = 10e-6;
    params.vec_init_simple();

    double gmax = 30e-3;
    Op_Gradient op_G(params, gmax);
    

    double smax = 200;
    Op_Slew op_S(params, smax);


    Op_Moment op_M0(params, 0, 2.0);
    Op_Moment op_M1(params, 1, 0.0);
    Op_Moment op_M2(params, 2, 0.0);

    op_M0.name = "M0";
    op_M1.name = "M1";  
    op_M2.name = "M2";
    
    Op_Identity op_I(params);

    // op_G.weight = 0.001;
    // op_S.weight = 0.001;
    // op_M0.weight = 1000.0;
    // op_M1.weight = 2000.0;
    // op_M2.weight = 2000.0;
    // op_I.obj_weight = -1000.0;


    op_G.weight = 100.0;
    op_S.weight = 0.00001;
    op_M0.weight = 100.0;
    op_M1.weight = 200.0;
    op_M2.weight = 200.0;
    op_I.obj_weight = -1.0;

    params.all_op.push_back(&op_G);
    params.all_op.push_back(&op_S);
    params.all_op.push_back(&op_M0);
    params.all_op.push_back(&op_M1);
    params.all_op.push_back(&op_M2);
    params.all_obj.push_back(&op_I);

    SolverGroptSDMM solver(params);
    solver.N_iter = 1000;
    solver.extra_debug = true;

    params.init();
    solver.solve();

    save_debug(params, solver);

    std::vector<double> vec(params.final_X.data(), params.final_X.data() + params.final_X.size());
    matplot::figure();
    matplot::plot(vec);

    std::vector<double> vec_slew(vec.size()-1, 0.0);
    for (int i = 0; i < vec_slew.size()-1; i++) {
        vec_slew[i] = (vec[i+1] - vec[i]) / params.dt;
    }

    matplot::figure();
    matplot::plot(vec_slew);


    matplot::show();

    spdlog::info("Finished main");

}


void diff_demo() 
{

    GroptParams params;
    params.diff_init_demo();

    double gmax = 30e-3;
    Op_Gradient op_G(params, gmax);
    
    double smax = 100;
    Op_Slew op_S(params, smax);

    Op_Moment op_M0(params, 0, 0.0);
    // Op_Moment op_M1(params, 1, 0.0);
    // Op_Moment op_M2(params, 2, 0.0);

    op_M0.name = "M0";
    // op_M1.name = "M1";  
    // op_M2.name = "M2";
    
    Op_BValue op_Bval(params);

    params.all_op.push_back(&op_G);
    params.all_op.push_back(&op_S);
    params.all_op.push_back(&op_M0);
    // params.all_op.push_back(&op_M1);
    // params.all_op.push_back(&op_M2);

    params.all_obj.push_back(&op_Bval);

    SolverGroptSDMM solver(params);
    solver.N_iter = 2000;
    solver.extra_debug = true;

    params.init();
    solver.solve();

    save_debug(params, solver);

    std::vector<double> vec(params.final_X.data(), params.final_X.data() + params.final_X.size());
    matplot::figure();
    matplot::plot(vec);

    std::vector<double> vec_slew(vec.size()-1, 0.0);
    for (int i = 0; i < vec_slew.size()-1; i++) {
        vec_slew[i] = (vec[i+1] - vec[i]) / params.dt;
    }

    matplot::figure();
    matplot::plot(vec_slew);


    matplot::show();

    spdlog::info("Finished main");

}


void diff_demo_c() 
{

    GroptParams params;
    params.diff_init_demo();

    double gmax = 30e-3;
    Op_Gradient op_G(params, gmax);
    
    double smax = 100;
    Op_Slew op_S(params, smax);

    Op_Moment op_M0(params, 0, 0.0);
    Op_Moment op_M1(params, 1, 0.0);
    Op_Moment op_M2(params, 2, 0.0);

    op_M0.name = "M0";
    op_M1.name = "M1";  
    op_M2.name = "M2";
    
    Op_BValue op_Bval(params);
    op_Bval.target = 300.0;
    op_Bval.tol0 = 50.0;

    params.all_op.push_back(&op_G);
    params.all_op.push_back(&op_S);
    params.all_op.push_back(&op_M0);
    params.all_op.push_back(&op_M1);
    params.all_op.push_back(&op_M2);
    params.all_op.push_back(&op_Bval);

    // Op_Identity op_I(params);
    // op_I.obj_weight = 1000.0;
    // params.all_obj.push_back(&op_I);

    SolverGroptSDMM solver(params);
    solver.N_iter = 10000;

    params.init();
    solver.solve();

    spdlog::info("bval = {:.1f}", op_Bval.get_bvalue(params.final_X));
    op_Bval.target = 310.0;
    
    // op_G.gmax = 29e-3;
    params.warm_start_prev();
    solver.solve();

    // op_G.gmax = 28.0e-3;
    // params.warm_start_prev();
    // solver.solve();

    spdlog::info("bval = {:.1f}", op_Bval.get_bvalue(params.final_X));

    save_debug(params, solver);

    std::vector<double> vec(params.final_X.data(), params.final_X.data() + params.final_X.size());
    matplot::figure();
    matplot::plot(vec);

    std::vector<double> vec_slew(vec.size()-1, 0.0);
    for (int i = 0; i < vec_slew.size()-1; i++) {
        vec_slew[i] = (vec[i+1] - vec[i]) / params.dt;
    }

    matplot::figure();
    matplot::plot(vec_slew);


    matplot::show();

    
    spdlog::info("Finished main");

}
*/

int main(int, char**){
    spdlog::set_level(spdlog::level::debug);

    // simple_test();
    // maximize_I();
    
    
}

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool


cdef extern from "gropt_params.hpp" namespace "Gropt":

    cdef cppclass GroptParams:

        GroptParams() except +

        int N
        int Naxis
        double dt
        int final_good
        int final_n_feval
        
        void vec_init_simple(int N, int Naxis, double first_val, double last_val)
        
        void diff_init(double _dt, double _TE, double _T_90, double _T_180, double _T_readout)

        void setvec_X0(int _N, int _Naxis, double *_X0, bool set_others)

        void set_ils_solver(string ils_method)

        void add_gmax(double gmax, bool rot_variant, double weight_mod)
        
        void add_smax(double smax, bool rot_variant, double weight_mod)
        
        void add_moment(double order, double target, 
                        double tol0, string units, int moment_axis, 
                        int start_idx0, int stop_idx0, int ref_idx0, double weight_mod)
        
        void add_SAFE(double stim_thresh,
                      double *tau1, double *tau2, double *tau3, 
                      double *a1, double *a2, double *a3,
                      double *stim_limit, double *g_scale,
                      int new_first_axis, bool demo_params, double weight_mod)

        void add_SAFE_vec(int N_vec, double *stim_thresh_vec,
                      double *tau1, double *tau2, double *tau3, 
                      double *a1, double *a2, double *a3,
                      double *stim_limit, double *g_scale,
                      int new_first_axis, bool demo_params, double weight_mod)
        
        void add_bvalue(double target, double tol, int start_idx0, int stop_idx0, double weight_mod)
        
        void add_TV(double tv_lam, double weight_mod)

        void add_obj_identity(double weight_mod)

        void prepare()
        void solve()
        void solve(int min_iter, 
                    int n_iter, 
                    double gamma_x, 
                    double ils_tol, 
                    int ils_max_iter, 
                    int ils_min_iter, 
                    double ils_sigma,
                    double ils_tik_lam) 

        void test_reduce_and_solve()

        void get_output(double **out, int &out_size)


cdef extern from "solver_groptsdmm.hpp" namespace "Gropt":
    cdef cppclass SolverGroptSDMM:

        SolverGroptSDMM() except +

        void set_general_params(int min_iter, int max_iter, int log_interval, double gamma_x, int max_feval)
        void set_ils_params(double ils_tol, int ils_max_iter, int ils_min_iter, double ils_sigma, double ils_tik_lam)
        void solve(GroptParams &gparams)
        void set_sdmm_params(int rw_interval, double rw_e_corr, double rw_eps, double rw_scalelim,
                             int grw_min_infeasible, int grw_interval, double grw_mod)

cdef extern from "gropt_utils.hpp" namespace "Gropt":
    
    void set_verbose(int level)
    void get_SAFE(int N, int Naxis, double dt, double *G_in, 
                  bool true_safe, int new_first_axis, bool demo_params,
                  double *tau1, double *tau2, double *tau3,
                  double *a1, double *a2, double *a3,
                  double *stim_limit, double *g_scale,
                  double **out, int &out_size)


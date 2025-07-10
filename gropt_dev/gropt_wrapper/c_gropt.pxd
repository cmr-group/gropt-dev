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

        
        void vec_init_simple(double first_val, double last_val)
        void diff_init(double _dt, double _TE, double _T_90, double _T_180, double _T_readout)

        void set_ils_solver(string ils_method)

        void add_gmax(double gmax)
        void add_smax(double smax)
        void add_moment(double order, double target)
        void add_SAFE(double stim_thresh,
                      double *tau1, double *tau2, double *tau3, 
                      double *a1, double *a2, double *a3,
                      double *stim_limit, double *g_scale,
                      int new_first_axis, bool demo_params)
        void add_bvalue(double target, double tol)
        void add_TV(double tv_lam, double weight_in)

        void add_obj_identity(double weight_mod)

        void init()
        void solve()
        void solve(int min_iter, 
                    int n_iter, 
                    double gamma_x, 
                    double ils_tol, 
                    int ils_max_iter, 
                    int ils_min_iter, 
                    double ils_sigma,
                    double ils_tik_lam) 

        void get_output(double **out, int &out_size)


cdef extern from "gropt_utils.hpp" namespace "Gropt":
    
    void set_verbose(int level)
    void get_SAFE(int N, int Naxis, double dt, double *G_in, 
                  bool true_safe, int new_first_axis, bool demo_params,
                  double *tau1, double *tau2, double *tau3,
                  double *a1, double *a2, double *a3,
                  double *stim_limit, double *g_scale,
                  double **out, int &out_size)


from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool


cdef extern from "gropt_params.hpp" namespace "Gropt":

    cdef cppclass GroptParams:

        GroptParams() except +

        int N
        int Naxis
        double dt

        void vec_init_simple()
        void diff_init_demo()

        void add_gmax(double gmax)
        void add_smax(double smax)
        void add_moment(double order, double target)
        void add_bvalue(double target, double tol)

        void add_obj_identity(double weight_mod)

        void init()
        void solve()
        void solve(int min_iter, 
                    int n_iter, 
                    double gamma_x, 
                    double ils_tol, 
                    int ils_max_iter, 
                    int ils_min_iter, 
                    double ils_sigma) 

        void get_output(double **out, int &out_size)

        
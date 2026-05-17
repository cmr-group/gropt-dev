#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/vector.h>
#include <nanobind/eigen/dense.h>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/callback_sink.h"

#include "gropt_params.hpp"
#include "solver_groptsdmm.hpp"
#include "solver_osqp.hpp"
#include "gropt_utils.hpp"
#include "equilibrate.hpp"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(gropt_wrapper, m) {
    m.doc() = "GrOpt: Gradient Optimization for MRI";
    m.attr("__build_date__") = __DATE__ " " __TIME__;


    // These allow you to get spdlogs into jupyter notebooks
    m.def("set_log_level", [](int level) {
      spdlog::set_level(static_cast<spdlog::level::level_enum>(level));
    }, "level"_a,
R"doc(Set the log level for the C++ gropt library.

Uses spdlog/Python logging convention: lower = more verbose.

Level mapping: 0=Trace, 1=Debug, 2=Info, 3=Warning, 4=Error, 5=Critical, 6=Off.

Parameters
----------
level : int
    Log level (0–6).)doc");

    m.def("set_log_callback", [](nb::callable cb) {
        auto sink = std::make_shared<spdlog::sinks::callback_sink_mt>(
            [cb](const spdlog::details::log_msg &msg) {
                std::string formatted(msg.payload.begin(), msg.payload.end());
                nb::gil_scoped_acquire gil;
                cb(static_cast<int>(msg.level), formatted);
            }
        );
        spdlog::default_logger()->sinks() = {sink};
    });

    
    //////////////////////////////////////////////////////////
    // -------------------------------------------------------
    // SolveResult
    // -------------------------------------------------------
    nb::class_<Gropt::SolveResult>(m, "SolveResult",
R"doc(Result from a GrOpt solve operation.

Attributes
----------
X : np.ndarray
    The optimized gradient waveform.
converged : bool
    Whether all constraints were satisfied.
n_iter : int
    Number of outer SDMM iterations.
n_feval : int
    Total number of inner linear solver iterations.)doc"
    )
        .def(nb::init<>())
        .def_rw("X", &Gropt::SolveResult::X)
        .def_rw("converged", &Gropt::SolveResult::converged)
        .def_rw("n_iter", &Gropt::SolveResult::n_iter)
        .def_rw("n_feval", &Gropt::SolveResult::n_feval)
        .def_rw("dt", &Gropt::SolveResult::dt)
        .def_rw("bvalue", &Gropt::SolveResult::bvalue)
        .def("__repr__", [](const Gropt::SolveResult &r) {
            return "SolveResult(converged=" + std::string(r.converged ? "True" : "False") +
                   ", n_iter=" + std::to_string(r.n_iter) +
                   ", n_feval=" + std::to_string(r.n_feval) +
                   ", X.size=" + std::to_string(r.X.size()) + ")";
        });
 
        
    //////////////////////////////////////////////////////////
    // -------------------------------------------------------
    // GroptParams
    // -------------------------------------------------------
    nb::class_<Gropt::GroptParams>(m, "GroptParams",
R"doc(Problem definition for GrOpt gradient optimization.

This class holds the waveform dimensions, constraints, and objectives
that define a gradient optimization problem.)doc"
    )
        .def(nb::init<>())

        // Properties: N, Naxis, dt (forwarded to pdata via reference members)
        .def_prop_rw("N",
            [](Gropt::GroptParams &self) { return self.N; },
            [](Gropt::GroptParams &self, int val) { self.N = val; },
            "Number of time points per axis.")
        .def_prop_rw("Naxis",
            [](Gropt::GroptParams &self) { return self.Naxis; },
            [](Gropt::GroptParams &self, int val) { self.Naxis = val; },
            "Number of gradient axes.")
        .def_prop_rw("dt",
            [](Gropt::GroptParams &self) { return self.dt; },
            [](Gropt::GroptParams &self, double val) { self.dt = val; },
            "Raster time in seconds.")

        // vec_init_simple
        .def("vec_init_simple", &Gropt::GroptParams::vec_init_simple,
            "N"_a = -1, "Naxis"_a = -1, "first_val"_a = 0.0, "last_val"_a = 0.0,
R"doc(Initialize the set_vec and inv_vec settings.

Parameters
----------
N : int, optional
    Number of points in a single axis. Negative values use existing value.
Naxis : int, optional
    Number of axes. Negative values use existing value.
first_val : float, optional
    Fixed value for the first point [mT/m].
last_val : float, optional
    Fixed value for the last point [mT/m].)doc"
        )

        // diff_init
        .def("diff_init", &Gropt::GroptParams::diff_init,
            "dt"_a = 400e-6, "TE"_a = 80e-3, "T_90"_a = 3e-3, "T_180"_a = 5e-3, "T_readout"_a = 16e-3,
R"doc(Initialize diffusion sequence parameters.

Parameters
----------
dt : float, optional
    Raster time in seconds.
TE : float, optional
    Echo time in seconds.
T_90 : float, optional
    Duration of excitation RF pulse in seconds (time from excitation,
    so half of full RF pulse duration).
T_180 : float, optional
    Duration of refocusing RF pulse in seconds.
T_readout : float, optional
    Time to TE of the readout in seconds.)doc"
        )

        // diff_init
        .def("diff_init_deadtime", &Gropt::GroptParams::diff_init_deadtime,
            "dt"_a = 400e-6, "TE"_a = 80e-3, "T_90"_a = 3e-3, "T_180"_a = 5e-3, "T_readout"_a = 16e-3,
R"doc(Initialize diffusion sequence parameters, but with forced deadtime to make "convnetional" waveforms.

Parameters
----------
dt : float, optional
    Raster time in seconds.
TE : float, optional
    Echo time in seconds.
T_90 : float, optional
    Duration of excitation RF pulse in seconds (time from excitation,
    so half of full RF pulse duration).
T_180 : float, optional
    Duration of refocusing RF pulse in seconds.
T_readout : float, optional
    Time to TE of the readout in seconds.)doc"
        )

        // diff_init_preencode
        .def("diff_init_preencode", &Gropt::GroptParams::diff_init_preencode,
            "dt"_a = 400e-6, "TE"_a = 80e-3, "T_90"_a = 3e-3, "T_180"_a = 5e-3, "T_readout"_a = 16e-3, "T_pre"_a = 0.0,
R"doc(Initialize diffusion sequence parameters with a pre-encoding period.

Parameters
----------
dt : float, optional
    Raster time in seconds.
TE : float, optional
    Echo time in seconds.
T_90 : float, optional
    Duration of excitation RF pulse in seconds (time from excitation,
    so half of full RF pulse duration).
T_180 : float, optional
    Duration of refocusing RF pulse in seconds.
T_readout : float, optional
    Time to TE of the readout in seconds.
T_pre : float, optional
    Duration of the pre-encoding period in seconds.

Returns
-------
int
    Number of pre-encoding time points (N_pre).)doc"
        )

        // setvec_X0 — accepts numpy array, infers N/Naxis from shape
        .def("setvec_X0", [](Gropt::GroptParams &self, nb::ndarray<double, nb::ndim<1>> X0, bool set_others) {
            Eigen::Map<const Eigen::VectorXd> x(X0.data(), X0.shape(0));
            self.setvec_X0(Eigen::VectorXd(x), 1, set_others);
        }, "X0"_a, "set_others"_a = true,
R"doc(Set the initial waveform guess (1D).

Parameters
----------
X0 : np.ndarray
    Initial guess for the waveform (warm start) [T/m].
set_others : bool, optional
    If True, update inv_vec, set_vals, and fixer with standard values.)doc"
        )

        .def("setvec_X0", [](Gropt::GroptParams &self, nb::ndarray<double, nb::ndim<2>> X0, bool set_others) {
            int Naxis = X0.shape(0);
            int N = X0.shape(1);
            // Flatten row-major 2D array to Eigen VectorXd
            Eigen::VectorXd flat(N * Naxis);
            const double *data = X0.data();
            for (int i = 0; i < N * Naxis; i++) {
                flat(i) = data[i];
            }
            self.setvec_X0(flat, Naxis, set_others);
        }, "X0"_a, "set_others"_a = true,
R"doc(Set the initial waveform guess (2D: Naxis x N).

Parameters
----------
X0 : np.ndarray
    Initial guess for the waveform (warm start) [T/m]. Shape (Naxis, N).
set_others : bool, optional
    If True, update inv_vec, set_vals, and fixer with standard values.)doc"
        )

        // setvec_set_vals — accepts numpy array, NaN = free, finite = fixed
        .def("setvec_set_vals", [](Gropt::GroptParams &self, nb::ndarray<double, nb::ndim<1>> set_vals) {
            Eigen::Map<const Eigen::VectorXd> v(set_vals.data(), set_vals.shape(0));
            self.setvec_set_vals(Eigen::VectorXd(v), 1);
        }, "set_vals"_a,
R"doc(Manually set the set_vals vector (1D), and update fixer from its NaN pattern.

NaN entries mark a point as free, finite entries lock the waveform to that value
at that index. The fixer mask is rebuilt automatically (1.0 = free, 0.0 = fixed),
and X0 is updated at the fixed positions to match.

Note: pdata.inv_vec is not touched — call vec_init_simple() (or one of the
diff_init variants) first if you need it allocated.

Parameters
----------
set_vals : np.ndarray
    1D array of length N. NaN = free, finite value = fixed.)doc"
        )

        .def("setvec_set_vals", [](Gropt::GroptParams &self, nb::ndarray<double, nb::ndim<2>> set_vals) {
            int Naxis = set_vals.shape(0);
            int N = set_vals.shape(1);
            Eigen::VectorXd flat(N * Naxis);
            const double *data = set_vals.data();
            for (int i = 0; i < N * Naxis; i++) {
                flat(i) = data[i];
            }
            self.setvec_set_vals(flat, Naxis);
        }, "set_vals"_a,
R"doc(Manually set the set_vals vector (2D: Naxis x N), and update fixer from its NaN pattern.

NaN entries mark a point as free, finite entries lock the waveform to that value
at that index. The fixer mask is rebuilt automatically, and X0 is updated at the
fixed positions to match.

Parameters
----------
set_vals : np.ndarray
    2D array of shape (Naxis, N). NaN = free, finite value = fixed.)doc"
        )

        // set_ils_solver
        .def("set_ils_solver", &Gropt::GroptParams::set_ils_solver,
            "ils_method"_a = "CG",
R"doc(Set the indirect solver method.

Parameters
----------
ils_method : str
    Solver method: 'CG', 'NLCG', or 'BiCGstabl' (case-sensitive).)doc"
        )

        // add_gmax
        .def("add_gmax", &Gropt::GroptParams::add_gmax,
            "gmax"_a = 0.03, "rot_variant"_a = true, "weight_mod"_a = 1.0,
R"doc(Add a maximum gradient amplitude constraint.

Parameters
----------
gmax : float, optional
    Maximum allowed gradient magnitude [T/m].
rot_variant : bool, optional
    If True, use rotationally invariant formulation.
weight_mod : float, optional
    Weighting factor for this constraint.)doc"
        )


        // add_concomitant
        .def("add_concomitant", &Gropt::GroptParams::add_concomitant,
            "start_idx"_a = 0, "rot_variant"_a = true, "weight_mod"_a = 1.0,
R"doc(Add a concomitant constraint.

Parameters
----------
start_idx : int, optional
    Starting index for the constraint (-1 = beginning).
rot_variant : bool, optional
    If True, use rotationally invariant formulation.
weight_mod : float, optional
    Weighting factor for this constraint.)doc"
        )

        // add_smax
        .def("add_smax", &Gropt::GroptParams::add_smax,
            "smax"_a = 80.0, "rot_variant"_a = true, "weight_mod"_a = 1.0,
R"doc(Add a maximum gradient slew rate constraint.

Parameters
----------
smax : float, optional
    Maximum allowed gradient slew rate [T/m/s].
rot_variant : bool, optional
    If True, use rotationally invariant formulation.
weight_mod : float, optional
    Weighting factor for this constraint.)doc"
        )

        // add_smax_vec
        .def("add_smax_vec", &Gropt::GroptParams::add_smax_vec,
            "smax_vec"_a, "rot_variant"_a = true, "weight_mod"_a = 1.0,
R"doc(Add a per-location slew rate constraint.

Parameters
----------
smax_vec : np.ndarray
    Slew-rate limit at each location [T/m/s]. Length (N-1) is broadcast
    across all axes; length Naxis*(N-1) sets a per-axis limit. Only
    supported with rot_variant=True.
rot_variant : bool, optional
    Must be True. A vector limit on the slew magnitude across axes is
    not supported.
weight_mod : float, optional
    Weighting factor for this constraint.)doc"
        )

        // add_moment
        .def("add_moment", &Gropt::GroptParams::add_moment,
            "order"_a = 0, "target"_a = 0.0, "tol"_a = 1e-6, "units"_a = "mT*ms/m",
            "axis"_a = 0, "start_idx"_a = -1, "stop_idx"_a = -1, "ref_idx"_a = 0,
            "weight_mod"_a = 1.0,
R"doc(Add a moment constraint.

Parameters
----------
order : int, optional
    Moment order.
target : float, optional
    Target moment value.
tol : float, optional
    Tolerance for satisfying the constraint.
units : str, optional
    Units: 'mT*ms/m', 'T*s/m', 'rad*s/m', or 's/m'.
axis : int, optional
    Axis for moment calculation.
start_idx : int, optional
    Starting index (-1 = beginning).
stop_idx : int, optional
    Stopping index (-1 = end).
ref_idx : int, optional
    Reference index (t=0 for moment calculations).
weight_mod : float, optional
    Weighting factor for this constraint.)doc"
        )

        // add_SAFE
        .def("add_SAFE", [](Gropt::GroptParams &self, double stim_thresh, int new_first_axis,
                            bool demo_params, nb::object safe_params_obj, double weight_mod) {
            if (!safe_params_obj.is_none()) {
                nb::dict sp = nb::cast<nb::dict>(safe_params_obj);
                Eigen::VectorXd tau1 = nb::cast<Eigen::VectorXd>(sp["tau1"]);
                Eigen::VectorXd tau2 = nb::cast<Eigen::VectorXd>(sp["tau2"]);
                Eigen::VectorXd tau3 = nb::cast<Eigen::VectorXd>(sp["tau3"]);
                Eigen::VectorXd a1 = nb::cast<Eigen::VectorXd>(sp["a1"]);
                Eigen::VectorXd a2 = nb::cast<Eigen::VectorXd>(sp["a2"]);
                Eigen::VectorXd a3 = nb::cast<Eigen::VectorXd>(sp["a3"]);
                Eigen::VectorXd stim_limit = nb::cast<Eigen::VectorXd>(sp["stim_limit"]);
                Eigen::VectorXd g_scale = nb::cast<Eigen::VectorXd>(sp["g_scale"]);
                self.add_SAFE(stim_thresh, tau1, tau2, tau3, a1, a2, a3, stim_limit, g_scale,
                              new_first_axis, weight_mod);
            } else {
                if (!demo_params) {
                    throw std::invalid_argument("If safe_params is None, demo_params must be True.");
                }
                self.add_SAFE(stim_thresh, new_first_axis, weight_mod);
            }
        }, "stim_thresh"_a = 1.0, "new_first_axis"_a = 0, "demo_params"_a = true,
           "safe_params"_a = nb::none(), "weight_mod"_a = 1.0,
R"doc(Add a SAFE (PNS) constraint.

Parameters
----------
stim_thresh : float, optional
    Stimulus threshold for the SAFE constraint.
new_first_axis : int, optional
    Swap the first axis of the SAFE parameters.
demo_params : bool, optional
    Whether to use demo parameters. Ignored if safe_params is provided.
safe_params : dict, optional
    Dictionary of SAFE parameters (see gropt.readasc).
weight_mod : float, optional
    Weighting factor for this constraint.)doc"
        )

        // add_SAFE_vec
        .def("add_SAFE_vec", [](Gropt::GroptParams &self, Eigen::VectorXd stim_thresh_vec,
                                int new_first_axis, bool demo_params,
                                nb::object safe_params_obj, double weight_mod) {
            if (!safe_params_obj.is_none()) {
                nb::dict sp = nb::cast<nb::dict>(safe_params_obj);
                Eigen::VectorXd tau1 = nb::cast<Eigen::VectorXd>(sp["tau1"]);
                Eigen::VectorXd tau2 = nb::cast<Eigen::VectorXd>(sp["tau2"]);
                Eigen::VectorXd tau3 = nb::cast<Eigen::VectorXd>(sp["tau3"]);
                Eigen::VectorXd a1 = nb::cast<Eigen::VectorXd>(sp["a1"]);
                Eigen::VectorXd a2 = nb::cast<Eigen::VectorXd>(sp["a2"]);
                Eigen::VectorXd a3 = nb::cast<Eigen::VectorXd>(sp["a3"]);
                Eigen::VectorXd stim_limit = nb::cast<Eigen::VectorXd>(sp["stim_limit"]);
                Eigen::VectorXd g_scale = nb::cast<Eigen::VectorXd>(sp["g_scale"]);
                self.add_SAFE_vec(stim_thresh_vec, tau1, tau2, tau3, a1, a2, a3,
                                  stim_limit, g_scale, new_first_axis, weight_mod);
            } else {
                if (!demo_params) {
                    throw std::invalid_argument("If safe_params is None, demo_params must be True.");
                }
                self.add_SAFE_vec(stim_thresh_vec, new_first_axis, weight_mod);
            }
        }, "stim_thresh_vec"_a, "new_first_axis"_a = 0, "demo_params"_a = true,
           "safe_params"_a = nb::none(), "weight_mod"_a = 1.0,
R"doc(Add a SAFE constraint with a vector stimulation limit.

Parameters
----------
stim_thresh_vec : np.ndarray
    Vector of stimulus thresholds.
new_first_axis : int, optional
    Swap the first axis of the SAFE parameters.
demo_params : bool, optional
    Whether to use demo parameters.
safe_params : dict, optional
    Dictionary of SAFE parameters (see gropt.readasc).
weight_mod : float, optional
    Weighting factor for this constraint.)doc"
        )

        // add_eddy
        .def("add_eddy", [](Gropt::GroptParams &self, nb::object lam_obj, double tol, double weight_mod) {
            Eigen::VectorXd lam;
            if (nb::isinstance<nb::float_>(lam_obj) || nb::isinstance<nb::int_>(lam_obj)) {
                lam.resize(1);
                lam(0) = nb::cast<double>(lam_obj);
            } else {
                lam = nb::cast<Eigen::VectorXd>(lam_obj);
            }
            self.add_eddy(lam, tol, weight_mod);
        }, "lam"_a, "tol"_a = .0001, "weight_mod"_a = 1.0,
R"doc(Add an eddy current constraint.

Parameters
----------
lam : float or np.ndarray
    Time constant(s) for eddy currents [seconds]. A single float is treated
    as a one-element array.
tol : float, optional
    Tolerance for the constraint.
weight_mod : float, optional
    Weighting factor for this constraint.)doc"
        )

        // add_bvalue
        .def("add_bvalue", [](Gropt::GroptParams &self, double target, double tol,
                              int start_idx0, int stop_idx0, double weight_mod,
                              nb::object mode_obj, double max_scale) {
            int mode_int;
            if (nb::isinstance<nb::str>(mode_obj)) {
                std::string mode_str = nb::cast<std::string>(mode_obj);
                if (mode_str == "setval") mode_int = 1;
                else if (mode_str == "minval") mode_int = 2;
                else if (mode_str == "minval_max") mode_int = 3;
                else throw std::invalid_argument("Invalid mode string. Must be 'setval', 'minval', or 'minval_max'.");
            } else {
                mode_int = nb::cast<int>(mode_obj);
                if (mode_int < 1 || mode_int > 3)
                    throw std::invalid_argument("Invalid mode integer. Must be 1, 2, or 3.");
            }
            self.add_bvalue(target, tol, start_idx0, stop_idx0, weight_mod, mode_int, max_scale);
        }, "target"_a = 100.0, "tol"_a = 1.0, "start_idx0"_a = -1, "stop_idx0"_a = -1,
           "weight_mod"_a = 1.0, "mode"_a = nb::int_(2), "max_scale"_a = 1.01,
R"doc(Add a b-value constraint.

Parameters
----------
target : float, optional
    Target b-value.
tol : float, optional
    Tolerance for the b-value constraint.
start_idx0 : int, optional
    Starting index (-1 = full waveform).
stop_idx0 : int, optional
    Stopping index (-1 = full waveform).
weight_mod : float, optional
    Weighting factor.
mode : int or str, optional
    'setval'/1 = set b-value, 'minval'/2 = minimum b-value (default),
    'minval_max'/3 = minimum b-value with scaling.
max_scale : float, optional
    Scale factor when mode=3.)doc"
        )

        // add_TV
        .def("add_TV", &Gropt::GroptParams::add_TV,
            "tv_lam"_a = 0.0, "weight_mod"_a = 1.0,
R"doc(Add total variation regularization.

Parameters
----------
tv_lam : float, optional
    Regularization strength (must be > 0 to have effect).
weight_mod : float, optional
    Weighting factor.)doc"
        )

        // add_obj_identity
        .def("add_obj_identity", &Gropt::GroptParams::add_obj_identity,
            "weight_mod"_a = 1.0,
R"doc(Add an identity (L2 norm) objective.

Parameters
----------
weight_mod : float, optional
    Weighting factor.)doc"
        )

        // prepare
        .def("prepare", &Gropt::GroptParams::prepare,
R"doc(Prepare the problem for solving.

Allocates vectors and sets up initial optimization variables.
Automatically called by solve() if not already done.)doc"
        )

        // print_op_details
        .def("print_op_details", &Gropt::GroptParams::print_op_details,
R"doc(Print details of all operators.

This function prints information about each operator in the problem, including their types and parameters.)doc"
        )

        
        // reset_op_weights
        .def("reset_op_weights", &Gropt::GroptParams::reset_op_weights,
R"doc(Reset all operator weights and spectral norms to 1.0.

Sets weight_mod, spec_norm, and spec_norm2 to 1.0 on every
operator in all_op and all_obj.)doc"
        );
    
    
    //////////////////////////////////////////////////////////
    // -------------------------------------------------------
    // Solver (base class)
    // -------------------------------------------------------
    nb::class_<Gropt::Solver>(m, "Solver")
        .def_rw("extra_debug", &Gropt::Solver::extra_debug,
            "If True, populate debug_solver with per-iteration history after solve().")
        .def("get_debug", [](Gropt::Solver &self) {
            nb::dict d;
            d["hist_X"]   = self.debug_solver.hist_X;
            d["hist_Ax"]  = self.debug_solver.hist_Ax;
            d["hist_z"]   = self.debug_solver.hist_z;
            d["hist_y"]   = self.debug_solver.hist_y;
            d["hist_Aty"] = self.debug_solver.hist_Aty;
            d["hist_weight"] = self.debug_solver.hist_weight;
            d["hist_gamma"] = self.debug_solver.hist_gamma;
            d["hist_gamma_x"] = self.debug_solver.hist_gamma_x;
            return d;
        },
R"doc(Return debug history as a dict of lists of numpy arrays.

Only populated when extra_debug=True before solve().

Returns
-------
dict with keys: hist_X, hist_Ax, hist_z, hist_y, hist_Aty
    Each value is a list of 1-D numpy arrays, one per logged iteration.)doc"
        )
        .def("set_general_params", &Gropt::Solver::set_general_params,
            "min_iter"_a = 1, "max_iter"_a = 2000, "log_interval"_a = 20,
            "gamma_x"_a = 1.6, "max_feval"_a = 12000, "extra_iters"_a = 0,
R"doc(Set general solver parameters.

Parameters
----------
min_iter : int, optional
    Minimum iterations.
max_iter : int, optional
    Maximum iterations.
log_interval : int, optional
    Logging interval (only visible with verbose logging).
gamma_x : float, optional
    Relaxation parameter for updates.
max_feval : int, optional
    Maximum total function evaluations.
extra_iters : int, optional
    Number of extra iterations to run after solution is found.)doc"
        )

        .def("set_ils_params", &Gropt::Solver::set_ils_params,
            "ils_tol"_a = 1e-3, "ils_max_iter"_a = 20, "ils_min_iter"_a = 2,
            "ils_sigma"_a = 1e-4, "ils_tik_lam"_a = 0.0,
R"doc(Set indirect linear solver parameters.

Parameters
----------
ils_tol : float, optional
    Relative tolerance for the inner solver.
ils_max_iter : int, optional
    Maximum ILS iterations per outer iteration.
ils_min_iter : int, optional
    Minimum ILS iterations.
ils_sigma : float, optional
    ADMM penalty parameter.
ils_tik_lam : float, optional
    Tikhonov regularization parameter.)doc"
        );


    //////////////////////////////////////////////////////////
    // -------------------------------------------------------
    // SolverGroptSDMM
    // -------------------------------------------------------
    nb::class_<Gropt::SolverGroptSDMM, Gropt::Solver>(m, "SolverGroptSDMM",
        "SDMM solver for GrOpt gradient optimization problems.")
        .def(nb::init<>())

        .def("set_sdmm_params", &Gropt::SolverGroptSDMM::set_sdmm_params,
            "rw_interval"_a = 8, "rw_e_corr"_a = 0.4, "rw_eps"_a = 1e-36,
            "rw_scalelim"_a = 1.5, "grw_min_infeasible"_a = 20,
            "grw_interval"_a = 20, "grw_mod"_a = 2.0,
R"doc(Set SDMM-specific parameters.

Parameters
----------
rw_interval : int, optional
    Interval for reweighting operations.
rw_e_corr : float, optional
    Error correction for reweighting.
rw_eps : float, optional
    Epsilon for numerical stability in reweighting.
rw_scalelim : float, optional
    Scale limit for reweighting.
grw_min_infeasible : int, optional
    Minimum infeasible iterations before adaptive reweighting.
grw_interval : int, optional
    Interval for adaptive reweighting checks.
grw_mod : float, optional
    Modification factor for adaptive reweighting.)doc"
        )

        .def("solve", [](Gropt::SolverGroptSDMM &self, Gropt::GroptParams &gparams) {
            Gropt::SolveResult result = self.solve(gparams);
            return result;
        }, "gparams"_a,
R"doc(Run the SDMM solver.

Parameters
----------
gparams : GroptParams
    The problem definition.

Returns
-------
SolveResult
    The optimization result containing the waveform and convergence info.)doc"
        );
    



    //////////////////////////////////////////////////////////
    // -------------------------------------------------------
    // SolverOSQP
    // -------------------------------------------------------
    nb::class_<Gropt::SolverOSQP, Gropt::Solver>(m, "SolverOSQP",
        "OSQP solver for GrOpt gradient optimization problems.")
        .def(nb::init<>())

        .def("solve", [](Gropt::SolverOSQP &self, Gropt::GroptParams &gparams) {
            Gropt::SolveResult result = self.solve(gparams);
            return result;
        }, "gparams"_a,
R"doc(Run the OSQP solver.

Parameters
----------
gparams : GroptParams
    The problem definition.

Returns
-------
SolveResult
    The optimization result containing the waveform and convergence info.)doc"
        );


    
    //////////////////////////////////////////////////////////
    // -------------------------------------------------------
    // Module-level functions
    // -------------------------------------------------------

    // solve() convenience function
    m.def("solve", [](Gropt::GroptParams &params,
                      int min_iter, int max_iter, int log_interval, double gamma_x, int max_feval, int extra_iters,
                      double ils_tol, int ils_max_iter, int ils_min_iter, double ils_sigma, double ils_tik_lam) {
        Gropt::SolverGroptSDMM solver;
        solver.set_general_params(min_iter, max_iter, log_interval, gamma_x, max_feval, extra_iters);
        solver.set_ils_params(ils_tol, ils_max_iter, ils_min_iter, ils_sigma, ils_tik_lam);
        return solver.solve(params);
    }, "params"_a,
       "min_iter"_a = 1, "max_iter"_a = 2000, "log_interval"_a = 20,
       "gamma_x"_a = 1.6, "max_feval"_a = 12000, "extra_iters"_a = 0,
       "ils_tol"_a = 1e-3, "ils_max_iter"_a = 20, "ils_min_iter"_a = 2,
       "ils_sigma"_a = 1e-4, "ils_tik_lam"_a = 0.0,
R"doc(Convenience function to solve a GrOpt problem.

Parameters
----------
params : GroptParams
    The problem definition.
min_iter, max_iter, log_interval, gamma_x, max_feval
    General solver parameters (see SolverGroptSDMM.set_general_params).
ils_tol, ils_max_iter, ils_min_iter, ils_sigma, ils_tik_lam
    ILS parameters (see SolverGroptSDMM.set_ils_params).

Returns
-------
SolveResult
    The optimization result.)doc"
    );

    // NormType enum
    nb::enum_<Gropt::NormType>(m, "NormType")
        .value("L2", Gropt::NormType::L2)
        .value("Inf", Gropt::NormType::Inf);

    // estimate_row_col_norms
    m.def("estimate_row_col_norms", [](Gropt::GroptParams &gparams, int n_reps, Gropt::NormType norm_type)
            -> std::pair<Eigen::VectorXd, Eigen::VectorXd> {
        Eigen::VectorXd row_norms, col_norms;
        Gropt::estimate_row_col_norms(gparams, n_reps, norm_type, row_norms, col_norms);
        return {row_norms, col_norms};
    }, "gparams"_a, "n_reps"_a = 10, "norm_type"_a = Gropt::NormType::Inf,
R"doc(Estimate row and column norms of the operator matrix.

Parameters
----------
gparams : GroptParams
    The problem definition (must have operators added).
n_reps : int, optional
    Number of random vector repetitions for the estimate.
norm_type : NormType, optional
    NormType.Inf (default) or NormType.L2.

Returns
-------
row_norms : np.ndarray
    Estimated norm for each constraint row.
col_norms : np.ndarray
    Estimated norm for each variable column.)doc"
    );

    // get_eq_vecs
    m.def("get_eq_vecs", [](Gropt::GroptParams &gparams)
            -> std::pair<Eigen::VectorXd, Eigen::VectorXd> {
        Eigen::VectorXd row_norms, col_norms;
        Gropt::get_eq_vecs(gparams, row_norms, col_norms);
        return {row_norms, col_norms};
    }, "gparams"_a,
R"doc(Get the current equilibration vectors from all operators.

Returns the accumulated eq_rows and eq_cols stored on each operator
after a call to equilibrate().

Parameters
----------
gparams : GroptParams
    The problem definition (must have operators prepared and equilibrated).

Returns
-------
row_norms : np.ndarray
    Accumulated row equilibration vector.
col_norms : np.ndarray
    Accumulated column equilibration vector.)doc"
    );

    // equilibrate
    m.def("equilibrate", &Gropt::equilibrate,
        "gparams"_a, "n_iter"_a = 5, "n_reps"_a = 10,
R"doc(Ruiz equilibration of the operator matrix.

Iteratively scales operator rows and columns so that the maximum
absolute value in each row and column approaches 1, improving
solver conditioning.

Parameters
----------
gparams : GroptParams
    The problem definition (operators will be modified in place).
n_iter : int, optional
    Number of equilibration iterations.
n_reps : int, optional
    Number of random vector repetitions per norm estimate.)doc"
    );

    // rescale_eq_vecs
    m.def("rescale_eq_vecs", &Gropt::rescale_eq_vecs,
        "gparams"_a, "row_scale"_a, "col_scale"_a,
R"doc(Rescale equilibration vectors by scalar factors.

Multiplies all operator eq_rows by row_scale and eq_cols by col_scale.

Parameters
----------
gparams : GroptParams
    The problem definition.
row_scale : float
    Scale factor applied to all row equilibration vectors.
col_scale : float
    Scale factor applied to all column equilibration vectors.)doc"
    );

    // estimate_spec_norm
    m.def("estimate_spec_norm", &Gropt::estimate_spec_norm,
        "gparams"_a, "n_iters"_a = 20,
R"doc(Estimate the spectral norm of the operator matrix via power iteration.

Parameters
----------
gparams : GroptParams
    The problem definition (must have operators prepared).
n_iters : int, optional
    Number of power iterations.

Returns
-------
float
    Estimated spectral norm.)doc"
    );

        // estimate_spec_norm
    m.def("estimate_individual_spec_norm", &Gropt::estimate_individual_spec_norm,
        "gparams"_a, "n_iters"_a = 20, "op_idx"_a = 0,
R"doc(Estimate the spectral norm of a single operator matrix via power iteration.

Parameters
----------
gparams : GroptParams
    The problem definition (must have operators prepared).
n_iters : int, optional
    Number of power iterations.
op_idx : int, optional
    Index of the operator for which to estimate the spectral norm.

Returns
-------
float
    Estimated spectral norm.)doc"
    );



    

    // get_SAFE
    m.def("get_SAFE", [](Eigen::VectorXd G, double dt, bool true_safe,
                         int new_first_axis, bool demo_params,
                         nb::object safe_params_obj) -> Eigen::VectorXd {
        int Naxis = 1;

        if (!safe_params_obj.is_none()) {
            nb::dict sp = nb::cast<nb::dict>(safe_params_obj);
            Eigen::VectorXd tau1 = nb::cast<Eigen::VectorXd>(sp["tau1"]);
            Eigen::VectorXd tau2 = nb::cast<Eigen::VectorXd>(sp["tau2"]);
            Eigen::VectorXd tau3 = nb::cast<Eigen::VectorXd>(sp["tau3"]);
            Eigen::VectorXd a1 = nb::cast<Eigen::VectorXd>(sp["a1"]);
            Eigen::VectorXd a2 = nb::cast<Eigen::VectorXd>(sp["a2"]);
            Eigen::VectorXd a3 = nb::cast<Eigen::VectorXd>(sp["a3"]);
            Eigen::VectorXd stim_limit = nb::cast<Eigen::VectorXd>(sp["stim_limit"]);
            Eigen::VectorXd g_scale = nb::cast<Eigen::VectorXd>(sp["g_scale"]);
            return Gropt::get_SAFE_eigen(G, Naxis, dt, true_safe, new_first_axis,
                                         tau1, tau2, tau3, a1, a2, a3, stim_limit, g_scale);
        } else {
            if (!demo_params) {
                throw std::invalid_argument("If safe_params is None, demo_params must be True.");
            }
            return Gropt::get_SAFE_eigen(G, Naxis, dt, true_safe, new_first_axis);
        }
    }, "G"_a, "dt"_a, "true_safe"_a = true, "new_first_axis"_a = 0,
       "demo_params"_a = true, "safe_params"_a = nb::none(),
R"doc(Compute the SAFE (PNS) response for a gradient waveform.

Parameters
----------
G : np.ndarray
    Gradient waveform.
dt : float
    Raster time in seconds.
true_safe : bool, optional
    Use true SAFE model.
new_first_axis : int, optional
    Swap the first axis of SAFE parameters.
demo_params : bool, optional
    Use demo parameters.
safe_params : dict, optional
    Dictionary of SAFE parameters (see gropt.readasc).

Returns
-------
np.ndarray
    SAFE response curve.)doc"
    );

    m.def("test_eigen_assertions", &Gropt::test_eigen_assertions,
        "test_type"_a,
R"doc(Test that assertions are running in Eigen. 

Each test is a different assertion that should be triggered in Eigen.  All tests are
expected to pass in normal execution.

Parameters
----------
test_type : int
    The type of assertion test to run, (1, 2, or 3).)doc"
    );

    
}

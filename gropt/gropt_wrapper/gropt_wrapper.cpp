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
#include "fft_tools.hpp"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(gropt_wrapper, m) {
    m.doc() = "GrOpt: Gradient Optimization for MRI";
    m.attr("__build_date__") = __DATE__ " " __TIME__;


    // Logging controls (gropt.setup_logging wraps these)
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
    },
R"doc(Route C++ log messages to a Python callable.

Replaces all sinks of the C++ default logger with one that calls fn(level, message),
where level is an int as in set_log_level and message is a str. Prefer
gropt.setup_logging(), which also releases the callback at exit.)doc");

    // Release the Python callback before shutdown; destroying it in C++ static teardown segfaults.
    m.def("clear_log_callback", []() {
        spdlog::default_logger()->sinks().clear();
    },
R"doc(Remove all C++ log sinks, releasing the Python callback; setup_logging registers this at exit.)doc");

    
    //////////////////////////////////////////////////////////
    // -------------------------------------------------------
    // SolveResult
    // -------------------------------------------------------
    nb::class_<Gropt::SolveResult>(m, "SolveResult",
R"doc(Result from a GrOpt solve operation.

Attributes
----------
X : np.ndarray
    The optimized gradient waveform [T/m], flat, length Naxis*N (axis-major).
converged : bool
    True if every constraint is feasible for the returned waveform.
n_iter : int
    Number of outer iterations.
n_feval : int
    Total number of inner linear-solver iterations.
dt : float
    Raster time of the waveform [s].
bvalue : float
    b-value of the returned waveform [s/mm^2] if the problem has a b-value term
    (constraint or objective), else 0.)doc"
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
        .def_rw("normalize_obj", &Gropt::GroptParams::normalize_obj,
            "If True, normalize the pull of each linearized objective (e.g. add_bvalue with "
            "as_objective=True) to a unit direction, so its magnitude is set by the objective weight "
            "(weight_mod) instead of growing with ||AᵀA x||. The weight then acts as a step-size "
            "knob. Default False.")

        .def_rw("eq_proj_solver", &Gropt::GroptParams::eq_proj_solver,
            "Linear solver for the equality projection used by project=True constraints. "
            "EqProjSolver.LDLT (default) is fast and precise but fails on singular or collinear rows. "
            "EqProjSolver.COD is a rank-revealing decomposition that handles near-collinear rows "
            "(e.g. several eddy time constants).")

        .def_rw("eq_proj_rcond", &Gropt::GroptParams::eq_proj_rcond,
            "Rank tolerance for EqProjSolver.COD: singular values below eq_proj_rcond times the largest "
            "are dropped; <= 0 uses Eigen's default. Ignored by LDLT. Default 1e-10.")

        .def_rw("safe_eps", &Gropt::GroptParams::safe_eps,
            "Smoothing [T/m/s] for the |.| in SAFE (PNS) operators: |v| -> sqrt(v^2 + eps^2), which "
            "makes their linearization continuous where the filtered slew crosses zero. 0 (default) = "
            "exact abs. Slightly conservative, since sqrt(v^2 + eps^2) >= |v|. Copied into each SAFE "
            "operator by add_SAFE/add_SAFE_vec, so set it first. A few percent of smax is a reasonable "
            "start (e.g. 1-5 for smax = 200).")

        // vec_init_simple
        .def("vec_init_simple", &Gropt::GroptParams::vec_init_simple,
            "N"_a = -1, "Naxis"_a = -1, "first_val"_a = 0.0, "last_val"_a = 0.0,
R"doc(Initialize a simple (non-diffusion) problem layout.

Sets inv_vec to +1, fixes the first and last points of every axis to first_val and
last_val (all other points free), and sets X0 to 0.01 at the free points.

Parameters
----------
N : int, optional
    Number of points in a single axis. Negative values use existing value.
Naxis : int, optional
    Number of axes. Negative values use existing value.
first_val : float, optional
    Fixed value for the first point [T/m].
last_val : float, optional
    Fixed value for the last point [T/m].)doc"
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

        // diff_init_deadtime
        .def("diff_init_deadtime", &Gropt::GroptParams::diff_init_deadtime,
            "dt"_a = 400e-6, "TE"_a = 80e-3, "T_90"_a = 3e-3, "T_180"_a = 5e-3, "T_readout"_a = 16e-3,
R"doc(Initialize diffusion sequence parameters with dead time for "conventional" waveforms.

Like diff_init, but the part of the pre-180 period that exceeds the post-180
period is fixed to zero (dead time), so both sides have the same free duration.

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
        .def("setvec_X0", [](Gropt::GroptParams &self, nb::ndarray<const double, nb::ndim<1>, nb::c_contig> X0,
                             bool set_others) {
            Eigen::Map<const Eigen::VectorXd> x(X0.data(), X0.shape(0));
            self.setvec_X0(Eigen::VectorXd(x), 1, set_others);
        }, "X0"_a, "set_others"_a = true,
R"doc(Set the initial waveform guess (1D).

Parameters
----------
X0 : np.ndarray
    Initial guess for the waveform (warm start) [T/m].
set_others : bool, optional
    If True, reset inv_vec to +1 and fix each axis's first and last points to
    their X0 values (all others free). Use False to keep a layout from diff_init.)doc"
        )

        .def("setvec_X0", [](Gropt::GroptParams &self, nb::ndarray<const double, nb::ndim<2>, nb::c_contig> X0,
                             bool set_others) {
            int Naxis = X0.shape(0);
            Eigen::Map<const Eigen::VectorXd> flat(X0.data(), X0.size());  // row-major (Naxis, N) -> axis-major
            self.setvec_X0(Eigen::VectorXd(flat), Naxis, set_others);
        }, "X0"_a, "set_others"_a = true,
R"doc(Set the initial waveform guess (2D: Naxis x N).

Parameters
----------
X0 : np.ndarray
    Initial guess for the waveform (warm start) [T/m]. Shape (Naxis, N).
set_others : bool, optional
    If True, reset inv_vec to +1 and fix each axis's first and last points to
    their X0 values (all others free). Use False to keep a layout from diff_init.)doc"
        )

        // setvec_set_vals — accepts numpy array, NaN = free, finite = fixed
        .def("setvec_set_vals", [](Gropt::GroptParams &self,
                                   nb::ndarray<const double, nb::ndim<1>, nb::c_contig> set_vals) {
            Eigen::Map<const Eigen::VectorXd> v(set_vals.data(), set_vals.shape(0));
            self.setvec_set_vals(Eigen::VectorXd(v), 1);
        }, "set_vals"_a,
R"doc(Manually set the set_vals vector (1D), and update fixer from its NaN pattern.

NaN entries mark a point as free, finite entries lock the waveform to that value
at that index. The fixer mask is rebuilt automatically (1.0 = free, 0.0 = fixed),
and X0 is updated at the fixed positions to match.

inv_vec is not changed; call vec_init_simple() (or one of the diff_init
variants) first if you need it allocated.

Parameters
----------
set_vals : np.ndarray
    1D array of length N. NaN = free, finite value = fixed.)doc"
        )

        .def("setvec_set_vals", [](Gropt::GroptParams &self,
                                   nb::ndarray<const double, nb::ndim<2>, nb::c_contig> set_vals) {
            int Naxis = set_vals.shape(0);
            Eigen::Map<const Eigen::VectorXd> flat(set_vals.data(), set_vals.size());
            self.setvec_set_vals(Eigen::VectorXd(flat), Naxis);
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

        // getvec_set_vals / getvec_fixer / getvec_inv_vec / getvec_X0 — flat 1D copies
        .def("getvec_set_vals", [](Gropt::GroptParams &self) { return Eigen::VectorXd(self.pdata.set_vals); },
R"doc(Return a copy of set_vals (flat, length Naxis*N): NaN = free, finite = fixed.)doc"
        )
        .def("getvec_fixer", [](Gropt::GroptParams &self) { return Eigen::VectorXd(self.pdata.fixer); },
R"doc(Return a copy of the fixer mask (flat, length Naxis*N): 1.0 = free, 0.0 = fixed.)doc"
        )
        .def("getvec_inv_vec", [](Gropt::GroptParams &self) { return Eigen::VectorXd(self.pdata.inv_vec); },
R"doc(Return a copy of inv_vec (flat, length Naxis*N): the +/- inversion pattern (sign flip at the 180).)doc"
        )
        .def("getvec_X0", [](Gropt::GroptParams &self) { return Eigen::VectorXd(self.pdata.X0); },
R"doc(Return a copy of the current initial-guess X0 (flat, length Naxis*N).)doc"
        )

        // set_ils_solver
        .def("set_ils_solver", &Gropt::GroptParams::set_ils_solver,
            "ils_method"_a = "CG",
R"doc(Set the inner (indirect) linear solver.

Only CG applies the equality projection used by project=True constraints; with
NLCG or BiCGstabl those constraints are enforced only by reproject_iterate.

Parameters
----------
ils_method : str
    Solver method: 'CG' (default), 'NLCG', or 'BiCGstabl' (case-sensitive).)doc"
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
    If True (default), limit each axis independently. If False, limit the
    gradient magnitude across axes, which is rotationally invariant.
weight_mod : float, optional
    Weighting factor for this constraint.)doc"
        )

        // add_gmax_vec
        .def("add_gmax_vec", &Gropt::GroptParams::add_gmax_vec,
            "gmax_vec"_a, "rot_variant"_a = true, "weight_mod"_a = 1.0,
R"doc(Add a per-location gradient amplitude constraint.

Parameters
----------
gmax_vec : np.ndarray
    Gradient amplitude limit at each location [T/m]. Length N is broadcast
    across all axes; length Naxis*N sets a per-axis limit. Only supported
    with rot_variant=True.
rot_variant : bool, optional
    Must be True. A vector limit on the gradient magnitude across axes is
    not supported.
weight_mod : float, optional
    Weighting factor for this constraint.)doc"
        )


        // add_concomitant
        .def("add_concomitant", &Gropt::GroptParams::add_concomitant,
            "start_idx"_a = 0, "rot_variant"_a = true, "weight_mod"_a = 1.0, "tol0"_a = 0.1,
            "target"_a = 1.0, "project"_a = false,
R"doc(Add a concomitant constraint.

Balances the gradient energy integral (sum of g^2 * dt) before and after the
180: constrains the ratio pos/neg to `target` (default 1 = balanced) within tol0.
This is a nonconvex quadratic constraint: by default a soft ADMM constraint, or
with project=True a relinearized equality projection.

The prox is the true Euclidean projection onto the feasible band: if the ratio
is inside the band it is left alone; otherwise the pre/post blocks are scaled by
minimal displacement onto the nearest boundary cone (reduces to the arithmetic
mean of the two block norms for exact equality).

Parameters
----------
start_idx : int, optional
    Index where the pre/post energy sums start (0 = beginning).
rot_variant : bool, optional
    Currently ignored; the energy balance is always computed across all axes.
weight_mod : float, optional
    Weighting factor for this constraint.
tol0 : float, optional
    Fractional tolerance on the pre/post energy ratio (feasible when
    |pos/neg - target| <= tol0). The prox aims for a cushioned (tighter) band so
    the solver reaches tol0 with margin. The dt in the energy integral cancels
    in this ratio for uniform dt.
target : float, optional
    Target pre/post energy ratio pos/neg (default 1.0 = balanced). The prox aims
    for pos/neg = target +/- tol; the projection enforces the linearized
    pos - target*neg = 0.
project : bool, optional
    If True, enforce the constraint via the equality null-space projection
    instead of as a soft ADMM constraint. The nonconvex energy balance is
    linearized at the current iterate each outer iteration (SQP) and projected
    jointly with any projected moments, so the prox is not used. weight_mod is
    then irrelevant and tol0 only sets the feasibility check.)doc"
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
    If True (default), limit each axis independently. If False, limit the
    slew magnitude across axes, which is rotationally invariant.
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
            "weight_mod"_a = 1.0, "project"_a = false, "absolute_tol"_a = false,
R"doc(Add a moment constraint.

Parameters
----------
order : int, optional
    Moment order (0, 1, 2, ...).
target : float, optional
    Target moment value, in `units`.
tol : float, optional
    Tolerance in M0 units (`units` at order 0). Order k uses tol * ||A_k|| / ||A_0||,
    the ratio of the moment row norms over the window, unless absolute_tol=True. With
    project=True the moment is projected onto the target and tol only sets the
    feasibility check.
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
    Weighting factor for this constraint.
project : bool, optional
    Enforce the moment via exact null-space projection instead of an ADMM penalty.
absolute_tol : bool, optional
    If True, tol is an absolute bound in this order's own units (e.g. for a
    nonzero M2 target). Default False (tol scaled from M0 as described above).)doc"
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
    Stimulation limit as a fraction of the SAFE threshold (1.0 = 100%).
new_first_axis : int, optional
    Use the SAFE parameters of this axis (0, 1, 2) for the first gradient axis
    (swapped with axis 0).
demo_params : bool, optional
    Use the built-in demo SAFE parameters. Must be True if safe_params is None;
    ignored if safe_params is provided.
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
    Per-sample stimulation limits as fractions of the SAFE threshold, length
    Naxis*N. Any other length falls back to 1.0 everywhere, with a warning.
new_first_axis : int, optional
    Use the SAFE parameters of this axis (0, 1, 2) for the first gradient axis
    (swapped with axis 0).
demo_params : bool, optional
    Use the built-in demo SAFE parameters. Must be True if safe_params is None;
    ignored if safe_params is provided.
safe_params : dict, optional
    Dictionary of SAFE parameters (see gropt.readasc).
weight_mod : float, optional
    Weighting factor for this constraint.)doc"
        )

        // add_eddy
        .def("add_eddy", [](Gropt::GroptParams &self, nb::object lam_obj, double tol, double weight_mod,
                            bool project) {
            Eigen::VectorXd lam;
            if (nb::isinstance<nb::float_>(lam_obj) || nb::isinstance<nb::int_>(lam_obj)) {
                lam.resize(1);
                lam(0) = nb::cast<double>(lam_obj);
            } else {
                lam = nb::cast<Eigen::VectorXd>(lam_obj);
            }
            self.add_eddy(lam, tol, weight_mod, project);
        }, "lam"_a, "tol"_a = .0001, "weight_mod"_a = 1.0, "project"_a = false,
R"doc(Add an eddy current constraint (residual eddy current at the end of the waveform).

The target is always 0. With project=False (default) it is a soft box constraint
|eddy| <= tol per time constant; with project=True the eddy term is held at 0 by
the null-space equality projection (like add_moment) inside the CG step instead
of an ADMM penalty.

Parameters
----------
lam : float or np.ndarray
    Time constant(s) for eddy currents [seconds]. A single float is treated
    as a one-element array. One row is added per (axis, time constant).
tol : float, optional
    Box half-width on the residual eddy term. With project=True it only sets the
    feasibility check.
weight_mod : float, optional
    Weighting factor for this constraint (constraint form only).
project : bool, optional
    If True, enforce eddy == 0 by null-space projection instead of the soft box
    constraint. Keep the number of time constants small, since many projected
    equalities (plus moments/concomitant) can over-constrain the free points.
    Close time constants give nearly collinear rows; use
    eq_proj_solver = EqProjSolver.COD in that case.)doc"
        )

        // add_bvalue
        .def("add_bvalue", [](Gropt::GroptParams &self, double target, double tol,
                              int start_idx0, int stop_idx0, double weight_mod,
                              nb::object mode_obj, double max_scale, bool as_objective, bool linearize) {
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
            self.add_bvalue(target, tol, start_idx0, stop_idx0, weight_mod, mode_int, max_scale, as_objective, linearize);
        }, "target"_a = 100.0, "tol"_a = 1.0, "start_idx0"_a = -1, "stop_idx0"_a = -1,
           "weight_mod"_a = 1.0, "mode"_a = nb::int_(2), "max_scale"_a = 1.01, "as_objective"_a = false,
           "linearize"_a = true,
R"doc(Add a b-value term (constraint by default, or a maximization objective).

Parameters
----------
target : float, optional
    Target b-value [s/mm^2]. Constraint only.
tol : float, optional
    Tolerance on the b-value [s/mm^2], used by mode='setval'.
start_idx0 : int, optional
    Starting index (-1 = full waveform).
stop_idx0 : int, optional
    Stopping index (-1 = full waveform).
weight_mod : float, optional
    Constraint: weighting factor. Objective: strength of the b-value pull
    (obj_weight = -weight_mod).
mode : int or str, optional
    Constraint only. 'setval'/1: b = target +/- tol. 'minval'/2 (default):
    b >= target. 'minval_max'/3: b >= target, and each prox also scales the
    waveform up by max_scale to keep pushing b higher.
max_scale : float, optional
    Per-iteration scale factor for mode='minval_max'.
as_objective : bool, optional
    If True, maximize the b-value as an objective instead of constraining it.
linearize : bool, optional
    Objective only. If True (default), the objective enters as a
    linearized gradient in the RHS (convex-concave/DCA), keeping the CG
    matrix positive-definite (recommended). If False, it enters as
    curvature in the LHS; this is experimental and can make the CG system
    indefinite, so watch for CG breakdown.)doc"
        )

        // add_diff_basin
        .def("add_diff_basin", &Gropt::GroptParams::add_diff_basin,
            "window_time"_a, "eps_factor"_a, "gmax"_a, "weight_mod"_a = 1.0, "same_sign"_a = false,
R"doc(Add a diffusion basin-orientation constraint (single 180 in the middle).

Requires the mean gradient in a short window on each side of the 180 to satisfy

    mean(g, pre window)  >= +eps
    mean(g, post window) <= -eps     (>= +eps with same_sign=True)

with eps = eps_factor * gmax. The 180 is located from the inv_vec sign flip, and
each window takes round(window_time / dt) free samples on its side, skipping the
fixed RF block. This selects which sign pattern b-value maximization converges
to and breaks its +/- symmetry. Keep eps_factor small so the constraint is
inactive at the optimum and does not bias the result. Enforced as a soft ADMM
constraint.

Parameters
----------
window_time : float
    Window width on each side of the 180 [s] (e.g. 1e-3).
eps_factor : float
    Minimum |mean window gradient| as a fraction of gmax (e.g. 0.07).
gmax : float
    Gradient limit [T/m], used to scale eps = eps_factor * gmax.
weight_mod : float, optional
    Weighting factor for this constraint.
same_sign : bool, optional
    False (default): opposite signs across the 180 (post <= -eps), the sign
    pattern of the global-max diffusion waveform. True: the same sign on both
    sides (post >= +eps).

Notes
-----
Assumes a single axis and a single 180. The pre window is always >= +eps; keep
any X0 seed orientation consistent with that.)doc"
        )

        // add_TV
        .def("add_TV", &Gropt::GroptParams::add_TV,
            "tv_lam"_a = 0.0, "weight_mod"_a = 1.0, "order"_a = 1,
R"doc(Add a total-variation (L1) penalty on a finite difference of the gradient.

Minimizes tv_lam * ||D g||_1, where D g is the physical slew [T/m/s] (order 1)
or jerk [T/m/s^2] (order 2). This is a penalty, not a constraint: it never
affects feasibility and is not used to select the returned iterate, so in a
problem with no other objective the result depends on min_iter.

Because the penalty is in physical units, useful tv_lam values are small and
depend on dt (e.g. ~1e-9 to 1e-6 for order 2 at dt = 400 us against a b-value
objective). Above that range the solution saturates at the minimum-TV waveform.

Parameters
----------
tv_lam : float, optional
    Penalty weight on ||D g||_1 in physical units; must be > 0 to have effect.
weight_mod : float, optional
    Initial ADMM weight (rho). Affects convergence, not the penalty strength.
order : int, optional
    1 = slew: sparse slew, blocky gradients; fights bang-bang ramps.
    2 = jerk: piecewise-constant slew; removes slew jitter without fighting
    the ramps.)doc"
        )

        // add_obj_identity
        .def("add_obj_identity", &Gropt::GroptParams::add_obj_identity,
            "weight_mod"_a = 1.0,
R"doc(Add an identity (L2 norm) objective that minimizes ||g||^2.

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

Prints each constraint and objective operator's type and parameters.)doc"
        )

        
        // reset_op_weights
        .def("reset_op_weights", &Gropt::GroptParams::reset_op_weights,
R"doc(Reset the spectral-norm estimates of all operators to 1.0.

Sets spec_norm and spec_norm2 to 1.0 on every operator in all_op and all_obj.
weight_mod is not changed.)doc"
        )

        // get_op_names
        .def("get_op_names", [](Gropt::GroptParams &self) {
            std::vector<std::string> names;
            for (auto &op : self.all_op) {
                names.push_back(op->name);
            }
            return names;
        },
R"doc(Return the constraint-operator names, in all_op order.

The order matches the per-operator index in the solver debug
histories (hist_weight, hist_gamma, hist_r_feas, hist_feas,
hist_con_pull_op), so the returned list can be used directly as plot
labels. hist_r_prim and hist_r_dual skip operators with project=True.)doc"
        )

        .def("get_op_keys", [](Gropt::GroptParams &self) {
            if (self.needs_prepare()) self.prepare(); // unique_name is assigned in prepare()
            std::vector<std::string> keys;
            for (auto &op : self.all_op) {
                keys.push_back(op->unique_name);
            }
            return keys;
        },
R"doc(Return each constraint operator's unique key ("<name>#<occurrence>"), in all_op order.

These are the keys warm-start snapshots use to match operators between solves.
Calls prepare() first if needed.)doc"
        );
    
    
    //////////////////////////////////////////////////////////
    // -------------------------------------------------------
    // Solver (base class)
    // -------------------------------------------------------
    nb::class_<Gropt::Solver>(m, "Solver")
        .def_rw("extra_debug", &Gropt::Solver::extra_debug,
            "If True, record per-iteration histories during solve(); read them with get_debug().")
        // --- General solver options (preferred over set_general_params) ---
        .def_rw("min_iter", &Gropt::Solver::min_iter,
            "Minimum outer iterations before feasible iterates are considered.")
        .def_rw("max_iter", &Gropt::Solver::max_iter, "Maximum outer iterations.")
        .def_rw("log_interval", &Gropt::Solver::log_interval, "Outer iterations between debug log prints.")
        .def_rw("gamma_x", &Gropt::Solver::gamma_x,
            "Over-relaxation factor for the outer primal update, X <- gamma_x*Xhat + (1-gamma_x)*X.")
        .def_rw("max_feval", &Gropt::Solver::max_feval,
            "Maximum total inner (CG) iterations over the whole solve.")
        .def_rw("obj_patience", &Gropt::Solver::obj_patience,
            "Objective problems: stop after this many feasible iters with no objective improvement.")
        .def_rw("obj_rtol", &Gropt::Solver::obj_rtol,
            "Relative objective-improvement threshold for the obj_patience plateau test.")
        // --- Inner linear-solver (ILS) options (preferred over set_ils_params) ---
        .def_rw("ils_tol", &Gropt::Solver::ils_tol,
            "Inner-solver tolerance, relative to the warm-start residual (~0.1 for inexact solves, "
            "smaller for near-exact).")
        .def_rw("ils_max_iter", &Gropt::Solver::ils_max_iter, "Max inner-solver iterations per outer step.")
        .def_rw("ils_min_iter", &Gropt::Solver::ils_min_iter, "Min inner-solver iterations per outer step.")
        .def_rw("ils_sigma", &Gropt::Solver::ils_sigma,
            "Proximal weight: adds sigma*I to the inner system, anchoring each inner solve to the "
            "current iterate.")
        .def_rw("ils_tik_lam", &Gropt::Solver::ils_tik_lam,
            "Tikhonov weight: adds tik_lam*I to the inner-system matrix, shrinking toward zero.")
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
            d["hist_r_prim"] = self.debug_solver.hist_r_prim;
            d["hist_r_dual"] = self.debug_solver.hist_r_dual;
            d["hist_r_feas"] = self.debug_solver.hist_r_feas;
            d["hist_feas"] = self.debug_solver.hist_feas;
            d["hist_all_feas"] = self.debug_solver.hist_all_feas;
            d["hist_bvalue"] = self.debug_solver.hist_bvalue;
            d["best_feasible_iter"] = self.debug_solver.best_feasible_iter;
            d["hist_cg_iter"] = self.debug_solver.hist_cg_iter;
            d["hist_cg_rnorm0"] = self.debug_solver.hist_cg_rnorm0;
            d["hist_cg_rnorm"] = self.debug_solver.hist_cg_rnorm;
            d["hist_cg_bnorm0"] = self.debug_solver.hist_cg_bnorm0;
            d["hist_obj_pull"] = self.debug_solver.hist_obj_pull;
            d["hist_con_pull"] = self.debug_solver.hist_con_pull;
            d["hist_con_pull_op"] = self.debug_solver.hist_con_pull_op;
            return d;
        },
R"doc(Return the debug histories of the last solve() as a dict.

Only populated when extra_debug=True before solve(). History lists have one
entry per outer iteration unless noted. Per-operator lists are ordered like
GroptParams.get_op_names(), except hist_r_prim and hist_r_dual, which skip
operators with project=True.

Returns
-------
dict with keys:
    hist_X, hist_Ax, hist_z, hist_y, hist_Aty
        list of 1-D numpy arrays.
    hist_weight, hist_gamma, hist_r_feas, hist_feas
        list of per-operator lists: ADMM weight, relaxation gamma, relative
        infeasibility, and binary feasible flag.
    hist_r_prim, hist_r_dual
        list of per-operator lists: Boyd ADMM primal/dual residuals.
    hist_all_feas
        list of ints: 1 if all operators were feasible that iteration, else 0.
    hist_gamma_x
        list of floats: the outer relaxation factor gamma_x.
    hist_bvalue
        b-value of the iterate (empty list if there is no b-value term).
    best_feasible_iter
        int: iteration of the returned iterate, or -1 if none was feasible
        (e.g. hist_bvalue[best_feasible_iter] is the returned b-value).
    hist_cg_iter, hist_cg_rnorm0, hist_cg_rnorm, hist_cg_bnorm0
        inner-solver diagnostics, one entry per inner solve (trust-region
        re-solves add entries): iterations taken, initial/final residual norm,
        and ||b||. hist_cg_iter[0] is a -1 placeholder for iteration 0.
    hist_obj_pull, hist_con_pull
        objective-vs-constraint balance: norm of the linearized-objective RHS
        pull and of the total constraint pull ||Σ Aᵀy||.
    hist_con_pull_op
        list of per-operator constraint pulls ||Aᵀy||.)doc"
        )
        .def("get_warmstart", [](Gropt::Solver &self) {
            Gropt::WarmStart w = self.get_warmstart();
            nb::dict d;
            d["active"] = w.active;
            d["N"] = w.N;
            d["Naxis"] = w.Naxis;
            d["dt"] = w.dt;
            d["X"] = Eigen::VectorXd(w.X);
            d["fixer"] = Eigen::VectorXd(w.fixer);
            nb::list ops;
            for (auto &o : w.ops) {
                nb::dict od;
                od["key"] = o.key;
                od["y"] = Eigen::VectorXd(o.y);
                od["weight"] = o.weight;
                od["gamma"] = o.gamma;
                od["spec_norm"] = o.spec_norm;
                od["blocks"] = o.blocks;
                ops.append(od);
            }
            d["ops"] = ops;
            return d;
        },
R"doc(Capture a warm-start snapshot from the last solve().

Returns the state at the returned (best-feasible) iterate, or at the final
iterate if none was feasible, as a plain, picklable dict. Pass it to
set_warmstart() before the next solve(). The consensus z is not stored; it is
regenerated as z = A x on load.

Returns
-------
dict with keys:
    active : bool   -- False if no snapshot is available (no solve has run).
    N, Naxis : int
    dt : float
    X, fixer : 1-D numpy arrays (length N*Naxis): primal and free/fixed mask.
    ops : list of dicts, one per operator (ordered like get_op_names()):
        key : str             -- operator key (see get_op_keys()), used for matching.
        y : 1-D numpy array   -- dual (Lagrange multiplier), in normalized Ax-space.
        weight, gamma : float -- ADMM weight and relaxation.
        spec_norm : float     -- operator normalization at capture, used to rescale y.
        blocks : list[int]    -- Ax-space partition, used to resize y.)doc"
        )
        .def("set_warmstart", [](Gropt::Solver &self, nb::dict d) {
            Gropt::WarmStart w;
            w.active = true;
            w.N = nb::cast<int>(d["N"]);
            w.Naxis = nb::cast<int>(d["Naxis"]);
            w.dt = nb::cast<double>(d["dt"]);
            w.X = nb::cast<Eigen::VectorXd>(d["X"]);
            w.fixer = nb::cast<Eigen::VectorXd>(d["fixer"]);
            for (nb::handle h : nb::cast<nb::list>(d["ops"])) {
                nb::dict od = nb::cast<nb::dict>(h);
                Gropt::OpWarmState st;
                st.key = nb::cast<std::string>(od["key"]);
                st.y = nb::cast<Eigen::VectorXd>(od["y"]);
                st.weight = nb::cast<double>(od["weight"]);
                st.gamma = nb::cast<double>(od["gamma"]);
                st.spec_norm = nb::cast<double>(od["spec_norm"]);
                st.blocks = nb::cast<std::vector<int>>(od["blocks"]);
                w.ops.push_back(st);
            }
            self.set_warmstart(w);
        }, "warmstart"_a,
R"doc(Load a warm-start snapshot (from get_warmstart()) for the next solve().

The snapshot applies to the next solve() only; later solves start cold (from
X0) unless this is called again.

The next solve() seeds the primal and each operator's dual, weight, and gamma
from the snapshot, resizing across a grid change: free runs of the waveform are
resampled, fixed runs come from the new set_vals, and each operator's dual is
resampled block by block. Operators are matched by key (see get_op_keys()); an
operator with no match (e.g. a newly added constraint) starts cold, so build the
new GroptParams in the same operator order. If Naxis or the number of free runs
per axis differs, the snapshot is ignored with a warning.

Parameters
----------
warmstart : dict
    A snapshot as returned by get_warmstart(). Requires keys N, Naxis, dt, X,
    fixer, and ops (each with key, y, weight, gamma, spec_norm, blocks).)doc"
        )
        .def("set_general_params", &Gropt::Solver::set_general_params,
            "min_iter"_a = 1, "max_iter"_a = 2000, "log_interval"_a = 20,
            "gamma_x"_a = 1.6, "max_feval"_a = 12000, "obj_patience"_a = 20,
R"doc(Set general solver parameters.

Deprecated: prefer setting the properties directly, e.g. `solver.max_iter = 5000`.

Parameters
----------
min_iter : int, optional
    Minimum outer iterations before feasible iterates are considered.
max_iter : int, optional
    Maximum outer iterations.
log_interval : int, optional
    Outer iterations between debug log prints.
gamma_x : float, optional
    Over-relaxation factor for the outer primal update.
max_feval : int, optional
    Maximum total inner (CG) iterations over the whole solve.
obj_patience : int, optional
    For objective problems, stop after this many feasible iterations with no
    objective improvement (returns the best feasible iterate). Ignored when
    there is no objective (returns on the first feasible iterate).)doc"
        )

        .def("set_ils_params", &Gropt::Solver::set_ils_params,
            "ils_tol"_a = 1e-3, "ils_max_iter"_a = 20, "ils_min_iter"_a = 2,
            "ils_sigma"_a = 1e-4, "ils_tik_lam"_a = 0.0,
R"doc(Set inner (indirect) linear solver parameters.

Deprecated: prefer setting the properties directly, e.g. `solver.ils_tol = 0.1`.

Parameters
----------
ils_tol : float, optional
    Inner-solver tolerance, relative to the warm-start residual.
ils_max_iter : int, optional
    Maximum inner-solver iterations per outer iteration.
ils_min_iter : int, optional
    Minimum inner-solver iterations per outer iteration.
ils_sigma : float, optional
    Proximal weight: adds sigma*I to the inner system.
ils_tik_lam : float, optional
    Tikhonov weight: adds tik_lam*I to the inner-system matrix.)doc"
        );


    //////////////////////////////////////////////////////////
    // -------------------------------------------------------
    // SolverGroptSDMM
    // -------------------------------------------------------
    nb::class_<Gropt::SolverGroptSDMM, Gropt::Solver>(m, "SolverGroptSDMM",
        "SDMM solver for GrOpt gradient optimization problems.")
        .def(nb::init<>())

        // --- SDMM reweighting options (preferred over set_sdmm_params) ---
        .def_rw("bb_enable", &Gropt::SolverGroptSDMM::bb_enable,
            "Per-operator BB (spectral) reweighting every rw_interval iterations (default True).")
        .def_rw("grw_enable", &Gropt::SolverGroptSDMM::grw_enable,
            "Global reweighting of the worst persistently infeasible operator (default True).")
        .def_rw("rw_interval", &Gropt::SolverGroptSDMM::rw_interval, "BB reweighting interval (iterations).")
        .def_rw("rw_e_corr", &Gropt::SolverGroptSDMM::rw_e_corr, "BB reweighting correlation threshold.")
        .def_rw("rw_eps", &Gropt::SolverGroptSDMM::rw_eps, "BB reweighting numerical-stability epsilon.")
        .def_rw("rw_scalelim", &Gropt::SolverGroptSDMM::rw_scalelim,
            "Max factor by which one BB update can raise or lower a weight.")
        .def_rw("grw_min_infeasible", &Gropt::SolverGroptSDMM::grw_min_infeasible,
            "grw: consecutive infeasible iterations before an operator can be bumped.")
        .def_rw("grw_interval", &Gropt::SolverGroptSDMM::grw_interval, "grw: reweighting interval (iterations).")
        .def_rw("grw_mod", &Gropt::SolverGroptSDMM::grw_mod, "grw: multiplicative weight-bump factor.")
        .def_rw("grw_balanced", &Gropt::SolverGroptSDMM::grw_balanced,
            "grw: if True, after bumping the worst constraint by grw_mod, divide every ADMM constraint "
            "weight by grw_mod^(1/K) (K = number of ADMM constraints) so their geometric mean stays "
            "fixed and the overall constraint scale relative to the objective is preserved. Default "
            "False (only the worst constraint is raised).")
        .def_rw("reproject_iterate", &Gropt::SolverGroptSDMM::reproject_iterate,
            "Re-project the over-relaxed iterate onto the equality (moment/eddy/concomitant) surface each "
            "outer iteration (default True). False lets the iterate drift off that surface between inner "
            "solves, which usually converges worse.")
        .def_rw("cutoff_freq", &Gropt::SolverGroptSDMM::cutoff_freq,
            "Low-frequency projection cutoff [Hz]; <= 0 disables. Each outer iteration the iterate is "
            "low-passed at cutoff_freq with a per-free-run DST-I to suppress high-frequency oscillation.")
        .def_rw("cutoff_iter", &Gropt::SolverGroptSDMM::cutoff_iter,
            "Outer iteration at which the cutoff_freq projection stops; < 0 = project on every iteration.")
        .def_rw("cutoff_trans", &Gropt::SolverGroptSDMM::cutoff_trans,
            "Raised-cosine roll-off width of the low-pass, as a fraction of the cutoff bin; 0 = brick wall "
            "(default). A wider roll-off suppresses more near-cutoff oscillation but also strips harmonics "
            "that keep plateaus flat; the net effect on plateau ripple depends on the waveform.")
        .def_rw("tr_enable", &Gropt::SolverGroptSDMM::tr_enable,
            "Trust-region step control (default False). After each inner solve, the tr_monitor signal "
            "checks whether the linearized model held over the step; if not, the step is re-solved from "
            "the current iterate with a larger proximal sigma. Keeps large objective steps from "
            "outrunning the nonlinear SAFE linearization.")
        .def_rw("tr_tol", &Gropt::SolverGroptSDMM::tr_tol,
            "Trust-region reject threshold; <= 0 uses the monitor default (0.2 for "
            "'linearization_error', 0.02 for 'feasibility', 0.5 for 'rel_step'). For "
            "'linearization_error' it is the max relative SAFE model error ||true - linear|| / ||true||.")
        .def_rw("tr_bump", &Gropt::SolverGroptSDMM::tr_bump,
            "Proximal-sigma multiplier applied on each rejected step (default 4).")
        .def_rw("tr_max_reject", &Gropt::SolverGroptSDMM::tr_max_reject,
            "Max re-solves per outer iteration before taking the most-damped step anyway (default 5).")
        .def_rw("tr_decay", &Gropt::SolverGroptSDMM::tr_decay,
            "Sigma relaxation factor toward ils_sigma on an accepted step (default 0.5).")
        .def_rw("tr_monitor", &Gropt::SolverGroptSDMM::tr_monitor,
            "Signal that drives the trust region: 'linearization_error' (default; relative error of the "
            "frozen SAFE linearization), 'feasibility' (reject if the SAFE violation grows past "
            "max(previous, tr_tol)), or 'rel_step' (||dx||/||x||). Only SAFE operators report the first "
            "two signals. Any other value disables the trust region with a warning.")
        .def_rw("obj_gate_enable", &Gropt::SolverGroptSDMM::obj_gate_enable,
            "Feasibility-gated objective (default False). Scales the objective pull by "
            "exp(-violation/obj_gate_scale), so the objective acts only near feasibility. Currently only "
            "SAFE (PNS/CNS) constraints report a violation; other constraints do not gate.")
        .def_rw("obj_gate_scale", &Gropt::SolverGroptSDMM::obj_gate_scale,
            "How sharply the objective gate opens as the constraint violation shrinks (in constraint "
            "units, e.g. ~0.05 of the SAFE limit). Smaller = stay gated closer to exact feasibility.")

        .def("set_sdmm_params", &Gropt::SolverGroptSDMM::set_sdmm_params,
            "rw_interval"_a = 8, "rw_e_corr"_a = 0.4, "rw_eps"_a = 1e-36,
            "rw_scalelim"_a = 1.5, "grw_min_infeasible"_a = 20,
            "grw_interval"_a = 20, "grw_mod"_a = 2.0,
R"doc(Set SDMM reweighting parameters.

Deprecated: prefer setting the properties directly, e.g. `solver.rw_interval = 16`.

Parameters
----------
rw_interval : int, optional
    BB reweighting interval (iterations).
rw_e_corr : float, optional
    BB reweighting correlation threshold.
rw_eps : float, optional
    BB reweighting numerical-stability epsilon.
rw_scalelim : float, optional
    Max factor by which one BB update can raise or lower a weight.
grw_min_infeasible : int, optional
    grw: consecutive infeasible iterations before an operator can be bumped.
grw_interval : int, optional
    grw: reweighting interval (iterations).
grw_mod : float, optional
    grw: multiplicative weight-bump factor.)doc"
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
                      int min_iter, int max_iter, int log_interval, double gamma_x, int max_feval, int obj_patience,
                      double ils_tol, int ils_max_iter, int ils_min_iter, double ils_sigma, double ils_tik_lam) {
        Gropt::SolverGroptSDMM solver;
        solver.set_general_params(min_iter, max_iter, log_interval, gamma_x, max_feval, obj_patience);
        solver.set_ils_params(ils_tol, ils_max_iter, ils_min_iter, ils_sigma, ils_tik_lam);
        return solver.solve(params);
    }, "params"_a,
       "min_iter"_a = 1, "max_iter"_a = 2000, "log_interval"_a = 20,
       "gamma_x"_a = 1.6, "max_feval"_a = 12000, "obj_patience"_a = 20,
       "ils_tol"_a = 1e-3, "ils_max_iter"_a = 20, "ils_min_iter"_a = 2,
       "ils_sigma"_a = 1e-4, "ils_tik_lam"_a = 0.0,
R"doc(Solve a GrOpt problem with a new SolverGroptSDMM.

Convenience wrapper that creates a SolverGroptSDMM with the given settings and
solves. For other options (reweighting, trust region, warm starts, debug
histories), create a SolverGroptSDMM and set its properties directly.

Parameters
----------
params : GroptParams
    The problem definition.
min_iter : int, optional
    Minimum outer iterations before feasible iterates are considered.
max_iter : int, optional
    Maximum outer iterations.
log_interval : int, optional
    Outer iterations between debug log prints.
gamma_x : float, optional
    Over-relaxation factor for the outer primal update.
max_feval : int, optional
    Maximum total inner (CG) iterations over the whole solve.
obj_patience : int, optional
    Objective problems: stop after this many feasible iterations with no
    objective improvement.
ils_tol : float, optional
    Inner-solver tolerance, relative to the warm-start residual.
ils_max_iter : int, optional
    Maximum inner-solver iterations per outer iteration.
ils_min_iter : int, optional
    Minimum inner-solver iterations per outer iteration.
ils_sigma : float, optional
    Proximal weight: adds sigma*I to the inner system.
ils_tik_lam : float, optional
    Tikhonov weight: adds tik_lam*I to the inner-system matrix.

Returns
-------
SolveResult
    The optimization result.)doc"
    );

    // NormType enum
    nb::enum_<Gropt::NormType>(m, "NormType")
        .value("L2", Gropt::NormType::L2)
        .value("Inf", Gropt::NormType::Inf);

    // EqProjSolver enum (GroptParams.eq_proj_solver)
    nb::enum_<Gropt::EqProjSolver>(m, "EqProjSolver")
        .value("LDLT", Gropt::EQ_LDLT)
        .value("COD", Gropt::EQ_COD);

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

    // estimate_individual_spec_norm
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
    Index of the operator in all_op (see GroptParams.get_op_names()).

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
R"doc(Compute the SAFE (PNS) response for a single-axis gradient waveform.

Parameters
----------
G : np.ndarray
    Gradient waveform [T/m], single axis.
dt : float
    Raster time in seconds.
true_safe : bool, optional
    Currently ignored.
new_first_axis : int, optional
    Use the SAFE parameters of this axis (0, 1, 2) for the waveform.
demo_params : bool, optional
    Use the built-in demo SAFE parameters. Must be True if safe_params is None.
safe_params : dict, optional
    Dictionary of SAFE parameters (see gropt.readasc).

Returns
-------
np.ndarray
    SAFE response at each time point, as a fraction of the stimulation limit
    (1.0 = at the limit).)doc"
    );

    // low_freq_project: exposes fft_tools LowFreqProjector (per-free-run DST-I low-pass)
    m.def("low_freq_project", [](Eigen::VectorXd x, double dt, double cutoff_hz,
                                 Eigen::VectorXd fixer, int Naxis, double trans_frac) -> Eigen::VectorXd {
        int N = static_cast<int>(x.size()) / Naxis;
        Gropt::LowFreqProjector proj;
        proj.setup(N, Naxis, dt, cutoff_hz, fixer, trans_frac);
        proj.project(x);
        return x;
    }, "x"_a, "dt"_a, "cutoff_hz"_a, "fixer"_a = Eigen::VectorXd(), "Naxis"_a = 1, "trans_frac"_a = 0.0,
R"doc(Low-pass a waveform with a per-free-run DST-I projection.

This is the filter SolverGroptSDMM applies when cutoff_freq > 0. Each maximal run
of free samples is band-limited independently with a DST-I, which assumes the
run is bounded by zeros; fixed samples are unchanged.

Parameters
----------
x : np.ndarray
    Waveform, length Naxis*N, axis-major.
dt : float
    Raster time [s].
cutoff_hz : float
    Cutoff frequency [Hz]; <= 0 returns x unchanged.
fixer : np.ndarray, optional
    Free mask (1 = free, 0 = fixed), length Naxis*N. Empty (default) or any
    other length treats every sample as free.
Naxis : int, optional
    Number of axes.
trans_frac : float, optional
    Raised-cosine roll-off width as a fraction of the cutoff bin; 0 (default) =
    brick wall. A wider roll-off suppresses more near-cutoff oscillation but
    also strips harmonics that keep plateaus flat.

Returns
-------
np.ndarray
    The projected copy of x.)doc"
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

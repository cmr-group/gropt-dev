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

    // Drop the callback sink so the captured Python callable is released NOW, while the interpreter is
    // still alive. The sink lives in the spdlog global default_logger; if that global is torn down during
    // C++ static destruction (after Py_Finalize) it destroys the captured nb::callable and decrefs a
    // Python object with no interpreter left -> SIGSEGV at exit. setup_logging registers this via atexit.
    m.def("clear_log_callback", []() {
        spdlog::default_logger()->sinks().clear();
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
        .def_rw("normalize_obj", &Gropt::GroptParams::normalize_obj,
            "If True, self-normalize the linearized objective direction so weight_mod sets a "
            "constant step magnitude (the pull no longer grows with ||AᵀA x||) — a forgiving rate "
            "knob for b-value maximization (endpoint is the constrained max regardless of value).")

        .def_rw("eq_proj_solver", &Gropt::GroptParams::eq_proj_solver,
            "Linear solver for the equality projection (project=True constraints). EqProjSolver.LDLT "
            "(default) is the original exact Gram solve -- fast/precise for well-conditioned sets "
            "(moments, concomitant) but blows up on singular/collinear rows. EqProjSolver.COD is a "
            "rank-revealing decomposition of Mhat that handles collinear rows (e.g. multiple eddy "
            "time-constants), keeps full precision (no Gram condition-squaring), and stays a true "
            "projection inside the CG. Switch to COD when LDLT can't handle your constraint set.")

        .def_rw("eq_proj_rcond", &Gropt::GroptParams::eq_proj_rcond,
            "Rank tolerance for eq_proj_solver = COD only (singular values below eq_proj_rcond * "
            "largest are dropped); ignored by LDLT. Larger drops more near-dependent rows; <= 0 uses "
            "Eigen's default threshold. Default 1e-10.")

        .def_rw("safe_eps", &Gropt::GroptParams::safe_eps,
            "Softabs smoothing (slew units, T/m/s) applied to every SAFE (PNS/CNS) op's |.|. 0 = hard "
            "abs (original). >0 replaces |v| with sqrt(v^2+eps^2) and the +-1 sign with a smooth "
            "v/sqrt(v^2+eps^2), so the frozen linearization changes continuously instead of flipping "
            "when a (filtered) slew crosses zero -- removes the near-zero sign churn that makes SAFE far "
            "more unstable than the linear eddy constraint. Slightly conservative. Must be set before "
            "add_SAFE. Try a few percent of smax (e.g. 1-5 for smax=200).")

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

        // getvec_set_vals / getvec_fixer / getvec_inv_vec / getvec_X0 — read-only COPIES (1D, flat)
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
            "target"_a = 1.0, "fix_gamma"_a = false, "gamma_fix"_a = 1.0, "project"_a = false,
            "as_objective"_a = false,
R"doc(Add a concomitant constraint.

Balances the gradient energy integral (sum of g^2 * dt) before and after the
180: constrains the ratio pos/neg to `target` (default 1 = balanced) within tol0.
This is a nonconvex quadratic constraint, so it is handled as a soft ADMM
constraint (it cannot use the linear moment null-space projection).

The prox is the true Euclidean projection onto the feasible band: if the ratio
is inside the band it is left alone; otherwise the pre/post blocks are scaled by
minimal displacement onto the nearest boundary cone (reduces to the arithmetic
mean of the two block norms for exact equality).

Parameters
----------
start_idx : int, optional
    Starting index for the constraint (-1 = beginning).
rot_variant : bool, optional
    If True, use rotationally invariant formulation.
weight_mod : float, optional
    Weighting factor for this constraint.
tol0 : float, optional
    Fractional tolerance on the pre/post energy ratio (feasible when
    |pos/neg - target| <= tol0). The prox aims for a cushioned (tighter) band so
    the solver reaches tol0 with margin. Note: the dt in the integral cancels
    in this ratio for uniform dt, so it does not change the numbers today --
    it keeps the energies physical and correct for non-uniform dt.
target : float, optional
    Target pre/post energy ratio pos/neg (default 1.0 = balanced). Applies in all
    modes: the box/prox aims for pos/neg = target +/- tol, and the projection /
    objective forms enforce the linearized pos - target*neg = 0.
fix_gamma : bool, optional
    If True, hold this operator's ADMM relaxation gamma fixed (the reweighter
    will not adapt it). The concomitant prox is nonconvex, so the reweighter's
    over-relaxation (gamma -> ~1.9) has no stability guarantee here; pin it to
    damp the resulting jitter.
gamma_fix : float, optional
    The fixed gamma used when fix_gamma is True. 1.0 = vanilla ADMM (no
    over-relaxation); < 1.0 = under-relaxed / extra damping.
    Ignored when project=True (no ADMM path).
project : bool, optional
    If True, enforce the constraint EXACTLY via the equality null-space
    projection instead of as a soft ADMM constraint. The nonconvex energy
    balance is linearized at the current iterate each outer iteration (SQP)
    and projected jointly with any projected moments, so it no longer uses
    the prox/consensus and the jitter from fighting the objective is removed.
    weight_mod/tol0/fix_gamma are then irrelevant (tol0 still sets the
    feasibility report band).
as_objective : bool, optional
    If True, drive the balance from the objective (all_obj) path via an
    augmented Lagrangian: c(x)=pos-neg is linearized each outer iteration and
    enforced with a scalar dual (method of multipliers) that injects a rank-1
    PSD curvature into the CG LHS and a pull into the RHS (weight_mod = the
    penalty rho). It contributes to the solve but is NOT the optimization
    target (does not affect best-feasible scoring) and is NOT feasibility-
    gated -- inspect pos/neg on the result. Mutually exclusive with project
    (project is ignored when as_objective=True). The three modes mirror
    b-value: default = soft ADMM constraint, project = exact projector,
    as_objective = augmented-Lagrangian objective.)doc"
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
            "weight_mod"_a = 1.0, "project"_a = false, "absolute_tol"_a = false,
R"doc(Add a moment constraint.

Parameters
----------
order : int, optional
    Moment order.
target : float, optional
    Target moment value.
tol : float, optional
    Order-0 (M0) feasibility tolerance. It is M0-anchored: higher-order moments scale their tolerance
    up by the row-norm ratio ||A_k|| / ||A_0|| (= (1e3*T_span)^k / sqrt(2k+1)), so this one number is
    order-consistent -- the same relative margin over each order's numerical floor -- in both projection
    and ADMM-box mode.
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
    Tolerance mode. False (default): M0-anchored -- `tol` is the order-0 tolerance and higher orders
    scale their bound by ||A_k||/||A_0||, so nulling is order-consistent. True: `tol` is an absolute
    bound in THIS order's physical units -- use with a nonzero higher-order `target` (e.g. a specified
    M2) when you want a fixed physical tolerance rather than a row-norm-scaled one.)doc"
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
|eddy| <= tol per time constant; with project=True the eddy current is driven to
EXACTLY 0 via the null-space equality projection (like add_moment), enforced
inside the CG step rather than as an ADMM penalty.

Parameters
----------
lam : float or np.ndarray
    Time constant(s) for eddy currents [seconds]. A single float is treated
    as a one-element array. One equality row is added per (axis, time constant).
tol : float, optional
    Box half-width for the constraint form (ignored when project=True).
weight_mod : float, optional
    Weighting factor for this constraint (constraint form only).
project : bool, optional
    If True, enforce eddy == 0 exactly by null-space projection instead of the
    soft box constraint. Use a small number of time constants -- many exact
    eddy equalities (plus moments/concomitant) can over-constrain the free DOFs.)doc"
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
    Scale factor when mode=3.
as_objective : bool, optional
    If True, add b-value as a maximization objective (all_obj,
    obj_weight=-1) instead of a constraint.
linearize : bool, optional
    Objective only. If True (default), the objective enters as a
    linearized gradient in the RHS (convex-concave/DCA), keeping the CG
    matrix positive-definite (recommended). If False, it enters as
    curvature in the LHS, which can make the system indefinite — for
    comparison only.)doc"
        )

        // add_diff_basin
        .def("add_diff_basin", &Gropt::GroptParams::add_diff_basin,
            "window_time"_a, "eps_factor"_a, "gmax"_a, "weight_mod"_a = 1.0, "same_sign"_a = false,
R"doc(Add a diffusion basin-orientation constraint (single 180 in the middle).

Forces the gradient to FLIP sign across the 180 -- the structure of the global-max
diffusion waveform -- and breaks the global +/- symmetry, so b-value maximization
commits to the correct basin in a SINGLE solve. Two one-sided constraints on the
mean gradient in a short window on each side of the 180 (located internally from
the inv_vec sign flip; the fixed 180 RF block is skipped):

    mean(g, pre window)  >= +eps
    mean(g, post window) <= -eps,    eps = eps_factor * gmax

Keep eps_factor small so the optimum (which has a strong flip) satisfies it with
slack -- the constraint is then inactive at the optimum and does not distort the
result; it only fences off the wrong basin during the descent. Soft ADMM
constraint, not an equality projection.

Parameters
----------
window_time : float
    Window width on each side of the 180 [s] (e.g. 1e-3).
eps_factor : float
    Min |mean window gradient| as a fraction of gmax (e.g. 0.07). Big enough to
    clear the ambiguous flat zone, small enough not to bias the lobe magnitude.
gmax : float
    Gradient limit [T/m], used to scale eps = eps_factor * gmax.
weight_mod : float, optional
    Weighting factor for this constraint.
same_sign : bool, optional
    False (default): force OPPOSITE signs across the 180 (pre >= +eps, post <= -eps)
    -- the flip / global-max basin. True: force the SAME sign (pre >= +eps,
    post >= +eps) -- the no-flip / local basin. (Negative eps does NOT do this; it
    only loosens the bound, so this is a flag.)

Notes
-----
Assumes a single axis / single 180 (the diffusion case). Default sign convention is
pre >= +eps, post <= -eps; keep any X0 seed orientation consistent with that.)doc"
        )

        // add_TV
        .def("add_TV", &Gropt::GroptParams::add_TV,
            "tv_lam"_a = 0.0, "weight_mod"_a = 1.0, "order"_a = 1,
R"doc(Add total-variation (L1) regularization on a finite difference of the gradient.

This is a PENALTY (prox of the L1 norm), not a budget constraint: it minimizes
tv_lam * ||difference||_1 rather than enforcing TV <= tol, and it does not gate
feasibility. tv_lam is a pure regularization WEIGHT on the physical slew
(order 1) / jerk (order 2): the prox divides by the operator's ADMM weight, so
the effective penalty is exactly tv_lam regardless of rho (the reweighter can
adapt rho for convergence without changing the result), and you do not have to
re-tune it when weight_mod changes. It is a weight, not a hard threshold -- use a
separate operator if you need a hard slew/jerk limit.

Parameters
----------
tv_lam : float, optional
    L1 regularization WEIGHT on the physical difference (must be > 0 to have
    effect). Larger = smoother.
weight_mod : float, optional
    Seeds this operator's initial ADMM weight (rho). Convergence/balance knob
    only -- it no longer changes the penalty strength (use tv_lam for that).
order : int, optional
    1 = TV of the gradient (||slew||_1): penalizes slew magnitude -> sparse
        slew / blocky gradients. Note this fights bang-bang ramps.
    2 = TV of the slew (||jerk||_1): penalizes CHANGES in slew -> piecewise-
        constant slew (clean bang-bang). This is the one that removes slew
        jitter without fighting the ramps; use a small tv_lam as a tie-breaker.)doc"
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
histories (hist_weight, hist_gamma, hist_r_prim, hist_r_dual),
so the returned list can be used directly as plot labels.)doc"
        );
    
    
    //////////////////////////////////////////////////////////
    // -------------------------------------------------------
    // Solver (base class)
    // -------------------------------------------------------
    nb::class_<Gropt::Solver>(m, "Solver")
        .def_rw("extra_debug", &Gropt::Solver::extra_debug,
            "If True, populate debug_solver with per-iteration history after solve().")
        // --- General solver options (preferred over set_general_params) ---
        .def_rw("min_iter", &Gropt::Solver::min_iter, "Minimum outer iterations before stopping.")
        .def_rw("max_iter", &Gropt::Solver::max_iter, "Maximum outer iterations.")
        .def_rw("log_interval", &Gropt::Solver::log_interval, "Iterations between debug log prints.")
        .def_rw("gamma_x", &Gropt::Solver::gamma_x, "Outer relaxation / over-relaxation factor.")
        .def_rw("max_feval", &Gropt::Solver::max_feval, "Maximum total inner (CG) iterations.")
        .def_rw("obj_patience", &Gropt::Solver::obj_patience,
            "Objective problems: stop after this many feasible iters with no objective improvement.")
        .def_rw("obj_rtol", &Gropt::Solver::obj_rtol,
            "Relative objective-improvement threshold for the obj_patience plateau test.")
        // --- Inner linear-solver (ILS) options (preferred over set_ils_params) ---
        .def_rw("ils_tol", &Gropt::Solver::ils_tol,
            "Inner-solver relative tolerance (warm-start-relative; ~0.1 for inexact, tight for near-exact).")
        .def_rw("ils_max_iter", &Gropt::Solver::ils_max_iter, "Max inner-solver iterations per outer step.")
        .def_rw("ils_min_iter", &Gropt::Solver::ils_min_iter, "Min inner-solver iterations per outer step.")
        .def_rw("ils_sigma", &Gropt::Solver::ils_sigma, "Proximal term added to the inner system (sigma I).")
        .def_rw("ils_tik_lam", &Gropt::Solver::ils_tik_lam, "Extra Tikhonov regularization for the inner system.")
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
            d["hist_cg_neg_curv"] = self.debug_solver.hist_cg_neg_curv;
            d["hist_cg_min_curv"] = self.debug_solver.hist_cg_min_curv;
            d["hist_cg_max_curv"] = self.debug_solver.hist_cg_max_curv;
            d["hist_obj_pull"] = self.debug_solver.hist_obj_pull;
            d["hist_con_pull"] = self.debug_solver.hist_con_pull;
            d["hist_con_pull_op"] = self.debug_solver.hist_con_pull_op;
            return d;
        },
R"doc(Return debug history as a dict of lists.

Only populated when extra_debug=True before solve(). Per-operator
lists are ordered to match GroptParams.get_op_names().

Returns
-------
dict with keys:
    hist_X, hist_Ax, hist_z, hist_y, hist_Aty
        list of 1-D numpy arrays, one per logged iteration.
    hist_weight, hist_gamma, hist_r_prim, hist_r_dual, hist_r_feas, hist_feas
        list (per iteration) of lists (per operator).
        hist_r_prim/hist_r_dual : Boyd ADMM primal/dual residuals.
        hist_r_feas/hist_feas   : relative infeasibility and binary feasible flag.
    hist_all_feas
        list of ints (per iteration): 1 if ALL operators were feasible that
        iteration, else 0. Plot against hist_bvalue to tell whether the
        objective is still climbing while feasibility flickers off (the
        returned best-feasible iterate is earlier) or has plateaued while
        feasible.
    hist_gamma_x
        list of floats, one per logged iteration.
    hist_bvalue
        achieved b-value per iteration (empty list if no b-value operator).
    best_feasible_iter
        int: the outer iteration whose iterate was actually returned (the
        best feasible one), or -1 if no feasible iterate was found. Use to
        mark the returned solution, e.g. hist_bvalue[best_feasible_iter].
    hist_cg_iter, hist_cg_rnorm0, hist_cg_rnorm, hist_cg_bnorm0
        inner linear-solver diagnostics, one entry per outer iteration:
        iterations taken, initial/final residual norm, and ||b||.
        (Generic to the ILS base class, so valid for CG/BiCGstabl/NLCG.)
    hist_cg_neg_curv
        CG only: 1 if a non-positive curvature direction (pAp<=0) was hit
        that outer iteration (indefinite system — watch this when running
        b-value as an objective with an unbalanced obj_weight).
    hist_cg_min_curv, hist_cg_max_curv
        CG only: min/max Rayleigh quotient (curvature) over the solve.
        hist_cg_min_curv -> 0 warns of near-singularity / blowup before
        hist_cg_neg_curv ever trips; max/min is a rough condition estimate.
        Only populated when extra_debug=True (the curvature dot is gated).
    hist_obj_pull, hist_con_pull
        objective-vs-constraint balance per iteration: ||g_obj|| (DCA
        objective RHS pull) and ||Σ Aᵀy|| (total constraint pull). The
        balance ratio is hist_obj_pull / hist_con_pull.
    hist_con_pull_op
        per-iteration list of per-operator ||Aᵀy|| (constraint pulls),
        ordered like get_op_names().)doc"
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
R"doc(Capture a warm-start snapshot from the last solve.

Returns the state at the returned (best-feasible) iterate as a plain, pickle-
friendly dict (e.g. to ship to loky workers). Pass it to set_warmstart() before
the next solve(). The consensus z is regenerated as z = A x on load, so it is
not stored.

Returns
-------
dict with keys:
    active : bool   -- False if no solve has run.
    N, Naxis : int
    dt : float
    X, fixer : 1-D numpy arrays (length N*Naxis): primal and source free/fixed mask.
    ops : list of dicts, one per operator (ordered like get_op_names()):
        key : str            -- operator unique_name, used for matching.
        y : 1-D numpy array  -- dual / Lagrange multiplier.
        weight, gamma : float
        blocks : list[int]   -- Ax-space partition (used to resize y).)doc"
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

One-shot: it applies to the next solve() only and is then cleared, so a later
solve() on the same solver is cold (uses setvec_X0 / X0) unless you call this
again. Chained sweeps still work since each step calls set_warmstart explicitly.

The next solve() seeds the primal, per-operator duals, and weights from it,
resizing across a grid change (free segments interpolate; each operator's Ax
blocks interpolate). Operators are matched by unique_name; an operator with no
match (e.g. a newly added constraint) starts cold. Build the new GroptParams in
the same operator order so the auto-assigned names line up.

Parameters
----------
warmstart : dict
    A snapshot as returned by get_warmstart().)doc"
        )
        .def("set_general_params", &Gropt::Solver::set_general_params,
            "min_iter"_a = 1, "max_iter"_a = 2000, "log_interval"_a = 20,
            "gamma_x"_a = 1.6, "max_feval"_a = 12000, "obj_patience"_a = 20,
R"doc([DEPRECATED -- prefer setting the properties directly, e.g. solver.max_iter = 5000] Set general solver parameters.

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
obj_patience : int, optional
    For objective problems, stop after this many feasible iterations with no
    objective improvement (returns the best feasible iterate). Ignored when
    there is no objective (returns on the first feasible iterate).)doc"
        )

        .def("set_ils_params", &Gropt::Solver::set_ils_params,
            "ils_tol"_a = 1e-3, "ils_max_iter"_a = 20, "ils_min_iter"_a = 2,
            "ils_sigma"_a = 1e-4, "ils_tik_lam"_a = 0.0,
R"doc([DEPRECATED -- prefer setting the properties directly, e.g. solver.ils_tol = 0.1] Set indirect linear solver parameters.

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

        // --- SDMM reweighting options (preferred over set_sdmm_params) ---
        .def_rw("rw_interval", &Gropt::SolverGroptSDMM::rw_interval, "BB reweighting interval (iterations).")
        .def_rw("rw_e_corr", &Gropt::SolverGroptSDMM::rw_e_corr, "BB reweighting correlation threshold.")
        .def_rw("rw_eps", &Gropt::SolverGroptSDMM::rw_eps, "BB reweighting numerical-stability epsilon.")
        .def_rw("rw_scalelim", &Gropt::SolverGroptSDMM::rw_scalelim, "BB reweighting per-step scale limit.")
        .def_rw("grw_min_infeasible", &Gropt::SolverGroptSDMM::grw_min_infeasible,
            "grw: minimum infeasible streak before adaptive reweighting kicks in.")
        .def_rw("grw_interval", &Gropt::SolverGroptSDMM::grw_interval, "grw: adaptive reweighting interval.")
        .def_rw("grw_mod", &Gropt::SolverGroptSDMM::grw_mod, "grw: multiplicative weight-bump factor.")
        .def_rw("grw_balanced", &Gropt::SolverGroptSDMM::grw_balanced,
            "grw: if True, REBALANCE instead of ratchet -- after bumping the worst constraint by grw_mod, "
            "divide every active constraint weight by grw_mod^(1/K) to hold their geometric mean fixed. "
            "Same emphasis shift on the worst, but the total constraint scale (vs the objective) is "
            "preserved, so the b-value pull isn't progressively drowned out.")
        .def_rw("reproject_iterate", &Gropt::SolverGroptSDMM::reproject_iterate,
            "Re-project the over-relaxed iterate onto the equality (moment/eddy/concomitant) surface every "
            "outer iteration (default True). Prevents the moment residual leaking under gamma_x != 1 + a "
            "loose CG. Set False to allow the old constraint-violating roaming (e.g. basin-crossing study).")
        .def_rw("cutoff_freq", &Gropt::SolverGroptSDMM::cutoff_freq,
            "Low-frequency projection cutoff [Hz]; <= 0 disables. Each outer iteration the iterate is "
            "projected onto frequencies <= cutoff_freq (per-axis DCT hard low-pass) to suppress "
            "high-frequency oscillation.")
        .def_rw("cutoff_iter", &Gropt::SolverGroptSDMM::cutoff_iter,
            "Outer iteration to STOP the cutoff_freq projection at; < 0 = project on every iteration.")
        .def_rw("cutoff_trans", &Gropt::SolverGroptSDMM::cutoff_trans,
            "Raised-cosine roll-off width of the low-pass, as a fraction of the cutoff bin "
            "(0 = brick wall, rings on plateaus; ~0.5 attenuates the near-cutoff oscillation band).")
        .def_rw("tr_enable", &Gropt::SolverGroptSDMM::tr_enable,
            "Trust-region step control (default False). After each inner CG, a StepMonitor checks whether "
            "the linearized model held over the step; if not, re-solve from X with the proximal sigma "
            "scaled up. Globalizes the nonlinear SAFE linearization so the objective's large steps can't "
            "outrun it (why SAFE blows up where a linear eddy/slew constraint holds).")
        .def_rw("tr_tol", &Gropt::SolverGroptSDMM::tr_tol,
            "Trust-region reject threshold. For tr_monitor='linearization_error' it is the max allowed "
            "relative SAFE model error ||true-linear||/||true|| (~0.1-0.3); <=0 keeps the monitor default.")
        .def_rw("tr_bump", &Gropt::SolverGroptSDMM::tr_bump,
            "Proximal-sigma multiplier applied on each rejected step (default 4).")
        .def_rw("tr_max_reject", &Gropt::SolverGroptSDMM::tr_max_reject,
            "Max re-solves per outer iteration before taking the most-damped step anyway (default 5).")
        .def_rw("tr_decay", &Gropt::SolverGroptSDMM::tr_decay,
            "Sigma relaxation factor toward ils_sigma on an accepted step (default 0.5).")
        .def_rw("tr_monitor", &Gropt::SolverGroptSDMM::tr_monitor,
            "Which divergence signal drives the trust region: 'linearization_error' (default, "
            "self-calibrating SAFE model fidelity), 'feasibility' (funnel), or 'rel_step' (||dx||/||x||).")
        .def_rw("obj_gate_enable", &Gropt::SolverGroptSDMM::obj_gate_enable,
            "Feasibility-gated objective (default False). Scales the objective pull by "
            "exp(-total_constraint_violation/obj_gate_scale): ~0 while infeasible (no cold overshoot past "
            "the constraints before they engage), ->1 when feasible (climb to max b). Fixes the fine-dt "
            "objective overshoot that no fixed bval_obj_weight can (too strong overshoots, too weak "
            "collapses).")
        .def_rw("obj_gate_scale", &Gropt::SolverGroptSDMM::obj_gate_scale,
            "How sharply the objective gate opens as the constraint violation shrinks (in constraint "
            "units, e.g. ~0.05 of the SAFE limit). Smaller = stay gated closer to exact feasibility.")

        .def("set_sdmm_params", &Gropt::SolverGroptSDMM::set_sdmm_params,
            "rw_interval"_a = 8, "rw_e_corr"_a = 0.4, "rw_eps"_a = 1e-36,
            "rw_scalelim"_a = 1.5, "grw_min_infeasible"_a = 20,
            "grw_interval"_a = 20, "grw_mod"_a = 2.0,
R"doc([DEPRECATED -- prefer setting the properties directly, e.g. solver.rw_interval = 16] Set SDMM-specific parameters.

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

    // resize_warmstart() -- run only the warm-start resizer, no solve (for inspecting the resize)
    m.def("resize_warmstart", [](Gropt::GroptParams &gp, nb::dict d) {
        // Parse the snapshot dict into a WarmStart (only the fields the resizer needs).
        Gropt::WarmStart w;
        w.active = true;
        w.N = nb::cast<int>(d["N"]);
        w.Naxis = nb::cast<int>(d["Naxis"]);
        w.X = nb::cast<Eigen::VectorXd>(d["X"]);
        w.fixer = nb::cast<Eigen::VectorXd>(d["fixer"]);
        for (nb::handle h : nb::cast<nb::list>(d["ops"])) {
            nb::dict od = nb::cast<nb::dict>(h);
            Gropt::OpWarmState st;
            st.key = nb::cast<std::string>(od["key"]);
            st.y = nb::cast<Eigen::VectorXd>(od["y"]);
            st.spec_norm = nb::cast<double>(od["spec_norm"]);
            st.blocks = nb::cast<std::vector<int>>(od["blocks"]);
            w.ops.push_back(st);
        }

        // The target problem must be prepared so operator sizes, masks, and unique_names exist.
        if (gp.op_prep_status != gp.N) gp.prepare();

        nb::dict out;
        bool compatible = (w.Naxis == gp.Naxis) &&
                          (Gropt::ws_free_run_count(w.fixer, w.Naxis) ==
                           Gropt::ws_free_run_count(gp.pdata.fixer, gp.Naxis));
        out["compatible"] = compatible;

        Eigen::VectorXd X = compatible ? Gropt::ws_resize_waveform(w.X, w.fixer, gp.pdata.fixer,
                                                                   gp.pdata.set_vals, gp.Naxis)
                                       : w.X;
        out["X"] = Eigen::VectorXd(X);

        nb::list ops;
        for (auto &op_ptr : gp.all_op) {
            Gropt::Operator *op = op_ptr.get();
            nb::dict od;
            od["key"] = op->unique_name;
            const Gropt::OpWarmState *st = w.find(op->unique_name);
            od["matched"] = (st != nullptr);
            if (st != nullptr) {
                od["y"] = Eigen::VectorXd(Gropt::ws_resize_dual(*st, op->Ax_block_lengths(), op->spec_norm));
            } else {
                od["y"] = Eigen::VectorXd(Eigen::VectorXd::Zero(op->Ax_size)); // unmatched -> cold
            }
            ops.append(od);
        }
        out["ops"] = ops;
        return out;
    }, "gparams"_a, "warmstart"_a,
R"doc(Resize a warm-start snapshot onto a target problem WITHOUT solving.

Runs only the resizer used inside solve(): segment-aware interpolation of the
primal X (free segments interpolate; fixed segments come from the target's
set_vals) and block-wise interpolation of each operator's dual y. Use it to
inspect the pre-optimization warm start; no equality projection or optimization
is applied.

Parameters
----------
gparams : GroptParams
    The TARGET problem (built at the new size). Prepared internally if needed.
warmstart : dict
    A snapshot from Solver.get_warmstart().

Returns
-------
dict with keys:
    compatible : bool   -- False if Naxis or free-segment topology differ (X then returned as-is).
    X : 1-D numpy array -- the resized primal on the target grid.
    ops : list of dicts (one per target operator, ordered like get_op_names()):
        key : str, matched : bool, y : 1-D numpy array (resized dual; zeros if unmatched).)doc"
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

    // low_freq_project: per-free-run DST-I low-pass projection -- exposes fft_tools.LowFreqProjector.
    // `fixer` is the binary free-mask (1=free, 0=fixed); pass an empty array to treat every sample free.
    m.def("low_freq_project", [](Eigen::VectorXd x, double dt, double cutoff_hz,
                                 Eigen::VectorXd fixer, int Naxis, double trans_frac) -> Eigen::VectorXd {
        int N = static_cast<int>(x.size()) / Naxis;
        Gropt::LowFreqProjector proj;
        proj.setup(N, Naxis, dt, cutoff_hz, fixer, trans_frac);
        proj.project(x);
        return x;
    }, "x"_a, "dt"_a, "cutoff_hz"_a, "fixer"_a = Eigen::VectorXd(), "Naxis"_a = 1, "trans_frac"_a = 0.0,
R"doc(Low-frequency projection of a waveform via per-free-run DST-I (fft_tools.LowFreqProjector).

x is length Naxis*N laid out axis-major. Each maximal run of free samples (fixer==1), bounded by fixed
zeros, is band-limited independently with a DST-I. The cutoff at cutoff_hz (per the dt grid) uses a
raised-cosine roll-off of fractional width trans_frac (0 = brick wall, the default). trans_frac>0 only
worsens plateau ripple (it strips the harmonics that flatten a plateau). Pass fixer empty to treat all
samples as free. Returns the projected copy.)doc"
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

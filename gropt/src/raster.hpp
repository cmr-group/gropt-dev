#ifndef GROPT_RASTER_H
#define GROPT_RASTER_H

/**
 * Exact conversion of a solved waveform to the hardware raster (e.g. 10 us). The scanner plays the
 * piecewise-linear function through the samples, so when dt_src = R * dt_tgt the target samples lie on the
 * same lines and the physical waveform is unchanged. Exact: gmax, slew, M0/M1 (zero endpoints), M2+ when
 * the lower moments are nulled or with pwl_moment_weight(), concomitant with Op_Concomitant::exact_quad,
 * SAFE with alpha_exact, b-value with Op_BValue::pwl_quad. With all four opted in, every constraint is
 * raster independent; with the defaults only the b-value moves, by ~gamma^2 (dt_src^2 - dt_tgt^2)/4 *
 * integral(g^2), ~0.1% at 400 us. Nothing is re-measured here (a second quadrature copy would drift):
 * certify with Operator::check().
 */

#include "Eigen/Dense"

namespace Gropt {

// Exact moment weight integral phi(t) t^k dt of the hat function at node time T, spacing h (one time unit):
// h, h*T for k < 2 (the rectangle rule), h*(T^2 + h^2/6), h*(T^3 + h^2 T/2), ... above. has_left/has_right:
// whether each neighbour is inside the window (an end node carries half a tent). Units (time unit)^(k+1).
double pwl_moment_weight(int k, double T, double h, bool has_left, bool has_right);

// R = dt_src / dt_tgt >= 1; throws std::invalid_argument unless it is an integer within rtol
int raster_ratio(double dt_src, double dt_tgt, double rtol = 1e-9);

// Linear resampling of X (axis-major, Naxis*N) onto an R-times finer raster: sample j -> j*R, natural
// length (N-1)*R + 1 per axis. A larger N_out zero-pads the end (e.g. when TE - T_readout is not a
// multiple of dt_src); a smaller one throws rather than drop gradient.
Eigen::VectorXd resample_waveform(const Eigen::VectorXd &X, int Naxis, int R, int N_out = -1);
Eigen::VectorXd resample_waveform(const Eigen::VectorXd &X, int Naxis, double dt_src, double dt_tgt,
                                  int N_out = -1);

// inv_vec for resample_waveform's result: each source sign held over its R target samples, so the flip
// stays at the same physical time; padding repeats the last sign.
Eigen::VectorXd resample_inv_vec(const Eigen::VectorXd &inv_vec, int Naxis, int R, int N_out = -1);

} // namespace Gropt

#endif

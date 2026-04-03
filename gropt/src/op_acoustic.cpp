#include "op_acoustic.hpp"
#include <cmath>
#include "spdlog/spdlog.h"

namespace Gropt {

Op_Acoustic::Op_Acoustic(const ProblemData &_pdata, std::vector<double> _freqs, std::vector<double> _bws, double _weight_mod,
                         double _bw_scale)
    : Operator(_pdata) {
    name = "Acoustic";
    freqs = _freqs;
    bws = _bws;
    weight_mod = _weight_mod;
    bw_scale = std::max(0.001, _bw_scale);
}

Op_Acoustic::Op_Acoustic(const Op_Acoustic& other)
    : Operator(other), freqs(other.freqs), bws(other.bws), bw_scale(other.bw_scale), H(other.H),
      N_pad(other.N_pad) {
    if (other.ffth) {
        ffth = std::make_unique<FFT_Helper>(other.N_pad);
    }
}

void Op_Acoustic::init() {
    // Dynamic padding size to dramatically speed up OSQP
    // Needs to be large enough for linear convolution without circular wrapping (+ some padding for high res representation)
    int n2 = pdata->N * 4;
    N_pad = 1024;
    while(N_pad < n2) N_pad *= 2;
    if (N_pad > 32768) N_pad = 32768;
    
    spdlog::trace("Op_Acoustic initialized with N_pad = {}", N_pad);

    Ax_size = pdata->Naxis * N_pad;

    H.setZero(N_pad);
    double df = 1.0 / (N_pad * pdata->dt);

    for (int k = 0; k < N_pad; k++) {
        double f = k * df;
        if (k > N_pad / 2) {
            f = (N_pad - k) * df;
        }

        double h_val = 0.0;
        for (size_t i = 0; i < freqs.size(); i++) {
            double dist = std::abs(f - freqs[i]);
            double bw = std::max(bws[i], 1e-6); // safeguard against zero division

            // Assume the target bandwidth represents the Full Width at Half Maximum (FWHM)
            // Scale the FWHM linearly by bw_scale
            double scaled_fwhm = bw * bw_scale;
            double sigma = scaled_fwhm / (2.0 * std::sqrt(2.0 * std::log(2.0)));
            double h_curr = std::exp(-0.5 * (dist * dist) / (sigma * sigma));

            if (h_curr > h_val) {
                h_val = h_curr;
            }
        }
        H(k) = h_val;
    }

    ffth = std::make_unique<FFT_Helper>(N_pad);

    spec_norm = 1.0;
    spec_norm2 = 1.0;

    if (do_init_weights) {
        obj_weight = 1e4; // Similar to b-value/slew constraints
        obj_weight *= weight_mod;
    }

    Operator::init();
}

void Op_Acoustic::forward(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    for (int i_ax = 0; i_ax < pdata->Naxis; i_ax++) {
        Eigen::VectorXd x_pad = Eigen::VectorXd::Zero(N_pad);
        x_pad.head(pdata->N) = X.segment(i_ax * pdata->N, pdata->N);
        Eigen::VectorXd out_pad(N_pad);
        ffth->fft_convolve(x_pad, out_pad, H, false, false);
        out.segment(i_ax * N_pad, N_pad) = out_pad;
    }
}

void Op_Acoustic::transpose(Eigen::VectorXd &X, Eigen::VectorXd &out) {
    // H is real and symmetric, so H^T = H
    for (int i_ax = 0; i_ax < pdata->Naxis; i_ax++) {
        Eigen::VectorXd x_pad = X.segment(i_ax * N_pad, N_pad);
        Eigen::VectorXd out_pad(N_pad);
        ffth->fft_convolve(x_pad, out_pad, H, false, true);
        out.segment(i_ax * pdata->N, pdata->N) = out_pad.head(pdata->N);
    }
}

void Op_Acoustic::prox(Eigen::VectorXd &X) {
    // Exact projection to A X = 0
    X.setZero();
}

void Op_Acoustic::check(Eigen::VectorXd &X) {
    // In strict sense, just check if norm of forbidden bands is near 0
    double err = X.norm();
    int feas = 1;
    if (err > 1e-3) {
        feas = 0;
    }
    hist_feas.push_back(feas);
}

} // namespace

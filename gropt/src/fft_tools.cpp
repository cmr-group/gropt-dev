#include "fft_tools.hpp"

#include <cmath>
#include <utility>

// Short waveforms: single-threaded (parallelism lives at a higher level) with a small plan cache.
#ifndef POCKETFFT_NO_MULTITHREADING
#define POCKETFFT_NO_MULTITHREADING
#endif
#ifndef POCKETFFT_CACHE_SIZE
#define POCKETFFT_CACHE_SIZE 8
#endif
#include "pocketfft_hdronly.hpp"

namespace Gropt {

static const double PI = std::acos(-1.0);

void LowFreqProjector::setup(int N, int Naxis, double dt, double cutoff_hz, const Eigen::VectorXd &fixer,
                             double trans_frac) {
    N_ = N;
    Naxis_ = Naxis;
    cutoff_hz_ = cutoff_hz;

    stride_ = {static_cast<std::ptrdiff_t>(sizeof(double))};
    axes_ = {0};
    shape_ = {0};
    coef_.assign(N_ > 0 ? N_ : 1, 0.0);

    runs_.clear();
    active_ = false;
    if (cutoff_hz <= 0.0 || dt <= 0.0 || N_ < 2) {
        return; // disabled / nothing to transform
    }

    const bool have_fixer = (fixer.size() == static_cast<Eigen::Index>(N_) * Naxis_);
    auto is_free = [&](int a, int local) -> bool {
        if (!have_fixer) return true;
        return fixer(static_cast<Eigen::Index>(a) * N_ + local) > 0.5;
    };

    // Split each axis into maximal free runs; each is low-passed independently.
    for (int a = 0; a < Naxis_; ++a) {
        int local = 0;
        while (local < N_) {
            if (!is_free(a, local)) { ++local; continue; }
            int start = local;
            while (local < N_ && is_free(a, local)) { ++local; }
            int M = local - start;

            int kc = static_cast<int>(std::floor(cutoff_hz * 2.0 * static_cast<double>(M + 1) * dt)) - 1;
            if (kc < 0) kc = 0;             // keep at least the fundamental half-sine
            if (kc > M - 1) kc = M - 1;     // cutoff at/above this run's Nyquist -> keep everything

            Run r;
            r.off = a * N_ + start;
            r.M = M;
            r.k_cut = kc;
            // pocketfft's DST-I (FFTW RODFT00) applied twice is 2*(M+1)*I; fct on both passes makes it
            // self-inverse.
            r.fct = 1.0 / std::sqrt(2.0 * static_cast<double>(M + 1));

            // Raised-cosine window (1 up to kpass, 0 at kc) to reduce brick-wall ringing on plateaus; a
            // wider taper rings less but rounds the passband more.
            r.win.assign(kc + 1, 1.0);
            int trans = static_cast<int>(std::lround((trans_frac > 0.0 ? trans_frac : 0.0) * (kc + 1)));
            if (trans > 0 && kc > 0) {
                int kpass = kc - trans;
                if (kpass < 0) kpass = 0;
                for (int k = kpass; k <= kc; ++k) {
                    double t = static_cast<double>(k - kpass) / static_cast<double>(kc - kpass);
                    r.win[k] = 0.5 * (1.0 + std::cos(PI * t)); // 1 at kpass -> 0 at kc
                }
            }

            runs_.push_back(std::move(r));
            if (kc < M - 1) active_ = true; // this run actually drops coefficients
        }
    }
    if (!active_) runs_.clear(); // no run has anything to remove -> full no-op
}

void LowFreqProjector::project(Eigen::VectorXd &x) const {
    if (!active_) {
        return;
    }
    for (const Run &r : runs_) {
        if (r.k_cut >= r.M - 1) {
            continue; // this run keeps every coefficient
        }
        shape_[0] = static_cast<std::size_t>(r.M);
        double *seg = x.data() + r.off;
        // forward DST-I seg -> coef_
        pocketfft::dst(shape_, stride_, stride_, axes_, 1, seg, coef_.data(), r.fct, false);
        // apply the window up to k_cut, zero above it
        for (int k = 0; k <= r.k_cut; ++k) {
            coef_[k] *= r.win[k];
        }
        for (int k = r.k_cut + 1; k < r.M; ++k) {
            coef_[k] = 0.0;
        }
        // inverse DST-I coef_ -> seg (DST-I is its own inverse with the same fct)
        pocketfft::dst(shape_, stride_, stride_, axes_, 1, coef_.data(), seg, r.fct, false);
    }
}

} // namespace Gropt

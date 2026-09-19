#ifndef FFT_TOOLS_H
#define FFT_TOOLS_H

/**
 * Frequency-domain tools backed by PocketFFT, which is included only in fft_tools.cpp.
 *
 * TODO: merge with FFT_Helper (fft_helper.hpp) from the GIRF work.
 */

#include <vector>

#include "Eigen/Dense"

namespace Gropt {

/**
 * Low-pass of a waveform's free samples via a DST-I per free run (x is axis-major, length Naxis*N).
 * Each maximal run of free samples is bounded by fixed zeros, so the sine basis matches that Dirichlet
 * boundary and the filtered run meets its fixed neighbours without a jump.
 * For a run of length M, bin k is at (k+1)/(2*(M+1)*dt) Hz; the last kept bin is
 * k_cut = floor(cutoff_hz*2*(M+1)*dt) - 1, clamped to [0, M-1].
 */
class LowFreqProjector {
  public:
    LowFreqProjector() = default;

    // fixer: free mask (1=free, 0=fixed), length N*Naxis, axis-major (GroptParams::pdata.fixer); on a size
    // mismatch every sample is free. The projector stays inactive (project() is a no-op) if cutoff_hz <= 0
    // or no run drops a coefficient. trans_frac: the window is 1 up to ~(1-trans_frac)*k_cut, then
    // cosine-tapers to 0 at k_cut; 0 (default) is a brick wall.
    void setup(int N, int Naxis, double dt, double cutoff_hz, const Eigen::VectorXd &fixer,
               double trans_frac = 0.0);

    // Low-pass x in place; no-op when inactive. Not thread-safe (shared mutable scratch).
    void project(Eigen::VectorXd &x) const;

    bool active() const { return active_; }
    int n_runs() const { return static_cast<int>(runs_.size()); }
    double cutoff_hz() const { return cutoff_hz_; }

  private:
    // One maximal run of free samples.
    struct Run {
        int off;    // absolute start index into x (axis*N + local start)
        int M;      // run length (number of free samples)
        int k_cut;  // last kept DST-I coefficient index (0 .. M-1)
        double fct; // 1/sqrt(2*(M+1)): makes DST-I its own inverse (pocketfft/FFTW RODFT00 convention)
        std::vector<double> win; // raised-cosine spectral weights for k = 0 .. k_cut (>k_cut is zeroed)
    };

    int N_ = 0;
    int Naxis_ = 1;
    double cutoff_hz_ = 0.0;
    bool active_ = false;

    std::vector<Run> runs_;

    // Plain-type equivalents of pocketfft::shape_t / stride_t, so this header needs no PocketFFT include.
    // shape_[0] is set to the run length before each transform.
    mutable std::vector<std::size_t> shape_;    // {M}
    std::vector<std::ptrdiff_t> stride_;         // {sizeof(double)}
    std::vector<std::size_t> axes_;             // {0}
    mutable std::vector<double> coef_;          // length-N DST scratch, reused per run
};

} // namespace Gropt

#endif

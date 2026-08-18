#ifndef FFT_TOOLS_H
#define FFT_TOOLS_H

/**
 * Frequency-domain tools (PocketFFT-backed). The PocketFFT header is included ONLY in fft_tools.cpp,
 * so nothing else in the build pulls it in -- the members below use the plain std::vector types that
 * pocketfft::shape_t / stride_t alias, so this header needs no PocketFFT include.
 * 
 * This is similar to the helper class from the GIRF work and needs to be merged together.
 */

#include <vector>

#include "Eigen/Dense"

namespace Gropt {

/**
 * Hard low-frequency projection via a per-free-run DST-I.
 *
 * A gradient waveform x (length Naxis*N, axis-major: [axis0(N), axis1(N), ...]) has FIXED samples
 * (fixer==0) and FREE samples (fixer==1). Each maximal run of free samples is bounded by fixed ZEROS on both
 * sides: a Dirichlet (zero-VALUE) boundary. The natural band-limited basis for such a segment is the DST
 * (sine) transform -- every sine basis function vanishes at the boundary, so a low-pass DST reconstruction
 * stays ~0 at the run edges and joins the fixed zeros with NO jump.
 *
 * Per free run of length M (free samples s..e, with x[s-1]=x[e+1]=0 fixed): DST-I -> zero every
 * coefficient above the cutoff bin -> DST-I (its own inverse). DST-I basis k has grid frequency
 *   f_k = (k+1) / (2*(M+1)*dt)   [Hz]   (k = 0 .. M-1)
 * so the last coefficient kept is k_cut = floor(cutoff_hz * 2*(M+1)*dt) - 1. Nyquist = 1/(2*dt).
 * Each run has its own M, hence its own bin resolution and k_cut, but the same physical cutoff_hz.
 */
class LowFreqProjector {
  public:
    LowFreqProjector() = default;

    // Configure for an (N, Naxis) problem on a dt grid, low-passing each free run at cutoff_hz. `fixer`
    // is the binary free-mask (1=free, 0=fixed), length N*Naxis, axis-major (GroptParams::pdata.fixer);
    // if its length doesn't match, every sample is treated as free (one run per axis). A cutoff <= 0, or
    // one that keeps every coefficient of every run, leaves the projector INACTIVE and project() a no-op.
    //
    // trans_frac sets the raised-cosine roll-off: the spectral window is 1 up to (1-trans_frac)*k_cut,
    // then cosine-tapers to 0 at k_cut (zero above). trans_frac=0 (default) is a brick wall.
    void setup(int N, int Naxis, double dt, double cutoff_hz, const Eigen::VectorXd &fixer,
               double trans_frac = 0.0);

    // Low-pass x in place (per free run). No-op when inactive. NOT thread-safe (shared scratch);
    // one projector per solver, single-threaded, as in the SDMM loop.
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

    // pocketfft::shape_t == std::vector<size_t>, stride_t == std::vector<ptrdiff_t>; kept as the plain
    // types so this header stays PocketFFT-free. shape_[0] is set to the run length before each transform.
    mutable std::vector<std::size_t> shape_;    // {M}
    std::vector<std::ptrdiff_t> stride_;         // {sizeof(double)}
    std::vector<std::size_t> axes_;             // {0}
    mutable std::vector<double> coef_;          // length-N DST scratch, reused per run
};

} // namespace Gropt

#endif

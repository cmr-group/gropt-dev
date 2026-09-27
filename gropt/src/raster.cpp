#include <cmath>
#include <stdexcept>
#include <string>

#include "raster.hpp"

namespace Gropt {

namespace {

// (1/h) * integral_0^h (h-u) (T+u)^k du -- the weight of the RIGHT half of the tent at T.
// Pass a negative h for the left half, negated: w_left(T,h,k) = -half_tent(k, T, -h).
double half_tent(int k, double T, double h) {
    double sum = 0.0;
    double binom = 1.0; // C(k,j), updated in place
    for (int j = 0; j <= k; j++) {
        sum += binom * std::pow(T, k - j) * std::pow(h, j + 1) / ((j + 1.0) * (j + 2.0));
        binom = binom * (k - j) / (j + 1.0);
    }
    return sum;
}

// Validates the sizes; N_in = source samples per axis, N_res = result samples per axis.
void resample_sizes(int in_size, int Naxis, int R, int N_out, int &N_in, int &N_res) {
    if (Naxis < 1) throw std::invalid_argument("resample: Naxis must be >= 1");
    if (R < 1) throw std::invalid_argument("resample: R must be >= 1");
    if (in_size < 1 || in_size % Naxis != 0) {
        throw std::invalid_argument("resample: input size must be a positive multiple of Naxis");
    }
    N_in = in_size / Naxis;
    const int N_nat = (N_in - 1) * R + 1;
    N_res = (N_out < 0) ? N_nat : N_out;
    if (N_res < N_nat) {
        throw std::invalid_argument("resample: N_out " + std::to_string(N_res) + " < natural length " +
                                    std::to_string(N_nat) + "; truncating would drop gradient samples");
    }
}

} // namespace

double pwl_moment_weight(int k, double T, double h, bool has_left, bool has_right) {
    if (k < 0) throw std::invalid_argument("pwl_moment_weight: k must be >= 0");
    double w = 0.0;
    if (has_right) w += half_tent(k, T, h);
    if (has_left) w -= half_tent(k, T, -h);
    return w;
}

int raster_ratio(double dt_src, double dt_tgt, double rtol) {
    if (!(dt_src > 0.0) || !(dt_tgt > 0.0)) {
        throw std::invalid_argument("raster_ratio: dt_src and dt_tgt must be > 0");
    }
    // dt_tgt > dt_src lands here too (R = 0, or R = 1 off by more than rtol)
    const double ratio = dt_src / dt_tgt;
    const int R = static_cast<int>(std::lround(ratio));
    if (R < 1 || std::fabs(ratio - R) > rtol * R) {
        throw std::invalid_argument("raster_ratio: dt_src is not an integer multiple of dt_tgt (ratio " +
                                    std::to_string(ratio) + ")");
    }
    return R;
}

Eigen::VectorXd resample_waveform(const Eigen::VectorXd &X, int Naxis, int R, int N_out) {
    int N_in, N_res;
    resample_sizes(static_cast<int>(X.size()), Naxis, R, N_out, N_in, N_res);

    Eigen::VectorXd out = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(Naxis) * N_res);
    for (int ax = 0; ax < Naxis; ax++) {
        const int si = ax * N_in;
        const int ti = ax * N_res;
        for (int j = 0; j < N_in - 1; j++) {
            const double a = X(si + j);
            const double b = X(si + j + 1);
            for (int m = 0; m < R; m++) { // m = 0 reproduces the source sample exactly
                const double f = static_cast<double>(m) / R;
                out(ti + j * R + m) = a + f * (b - a);
            }
        }
        out(ti + (N_in - 1) * R) = X(si + N_in - 1);
        // samples past the natural length stay zero (the waveform already ends at zero)
    }
    return out;
}

Eigen::VectorXd resample_waveform(const Eigen::VectorXd &X, int Naxis, double dt_src, double dt_tgt, int N_out) {
    return resample_waveform(X, Naxis, raster_ratio(dt_src, dt_tgt), N_out);
}

Eigen::VectorXd resample_inv_vec(const Eigen::VectorXd &inv_vec, int Naxis, int R, int N_out) {
    int N_in, N_res;
    resample_sizes(static_cast<int>(inv_vec.size()), Naxis, R, N_out, N_in, N_res);

    Eigen::VectorXd out(static_cast<Eigen::Index>(Naxis) * N_res); // every sample is written below
    for (int ax = 0; ax < Naxis; ax++) {
        const int si = ax * N_in;
        const int ti = ax * N_res;
        for (int j = 0; j < N_in; j++) { // hold each source sign over its R target samples; the last pads
            const double s = inv_vec(si + j);
            const int lo = j * R;
            const int hi = (j == N_in - 1) ? N_res : (j + 1) * R;
            for (int i = lo; i < hi; i++) out(ti + i) = s;
        }
    }
    return out;
}

} // namespace Gropt

#ifndef GROPT_WARMSTART_H
#define GROPT_WARMSTART_H

// Warm-starting one solve from another (e.g. a TE sweep, or adding a constraint). Carries the primal X
// and each operator's dual y, weight and gamma, matched by Operator::unique_name ("<name>#<occurrence>").
// z is regenerated as A*X on load; unmatched operators start cold.

#include "Eigen/Dense"
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace Gropt {

// Snapshot of one operator's ADMM state.
struct OpWarmState {
    std::string key;         // Operator::unique_name
    Eigen::VectorXd y;       // dual, in the operator's normalized Ax-space
    double weight = 1.0;     // ADMM penalty rho_i
    double gamma = 1.5;      // relaxation gamma_i
    double spec_norm = 1.0;  // operator spec_norm at capture, used to rescale y on load
    std::vector<int> blocks; // partition of y at capture (e.g. SAFE: n_terms*Naxis blocks of N)
};

// Full warm-start snapshot, decoupled from the operators/params that produced it.
struct WarmStart {
    bool active = false;
    int N = 0;
    int Naxis = 0;
    double dt = 0.0;
    Eigen::VectorXd X;     // primal at capture (length N*Naxis)
    Eigen::VectorXd fixer; // source free(!=0)/fixed(==0) mask (length N*Naxis)
    std::vector<OpWarmState> ops;

    const OpWarmState *find(const std::string &key) const {
        for (const auto &o : ops) {
            if (o.key == key) return &o;
        }
        return nullptr;
    }
};

// ---- resize primitives ---------------------------------------------------- //

// Linearly resample a 1-D vector to n_new samples, preserving the endpoints.
inline Eigen::VectorXd ws_resample(const Eigen::VectorXd &v, int n_new) {
    const int n = static_cast<int>(v.size());
    if (n_new <= 0) return Eigen::VectorXd();
    if (n_new == n) return v;
    if (n <= 1) return Eigen::VectorXd::Constant(n_new, n == 1 ? v(0) : 0.0);
    if (n_new == 1) { // the endpoint map below divides by n_new - 1; take the midpoint instead
        const double t = 0.5 * (n - 1);
        const int lo = static_cast<int>(t);
        const int hi = std::min(lo + 1, n - 1);
        return Eigen::VectorXd::Constant(1, (1.0 - (t - lo)) * v(lo) + (t - lo) * v(hi));
    }
    Eigen::VectorXd out(n_new);
    for (int i = 0; i < n_new; i++) {
        const double t = static_cast<double>(i) * (n - 1) / (n_new - 1); // [0,n_new-1] -> [0,n-1]
        const int lo = static_cast<int>(std::floor(t));
        const int hi = std::min(lo + 1, n - 1);
        const double f = t - lo;
        out(i) = (1.0 - f) * v(lo) + f * v(hi);
    }
    return out;
}

// Resample each block of v (partition blocks_src) to its blocks_tgt length independently, so stacked
// blocks never bleed into each other. Returns an empty vector if the partitions are incompatible.
inline Eigen::VectorXd ws_resize_blocks(const Eigen::VectorXd &v, const std::vector<int> &blocks_src,
                                        const std::vector<int> &blocks_tgt) {
    if (blocks_src.size() != blocks_tgt.size()) return Eigen::VectorXd();
    int src_total = 0, tgt_total = 0;
    for (int b : blocks_src) src_total += b;
    for (int b : blocks_tgt) tgt_total += b;
    if (src_total != v.size()) return Eigen::VectorXd();
    Eigen::VectorXd out(tgt_total);
    int si = 0, ti = 0;
    for (size_t k = 0; k < blocks_src.size(); k++) {
        out.segment(ti, blocks_tgt[k]) = ws_resample(v.segment(si, blocks_src[k]), blocks_tgt[k]);
        si += blocks_src[k];
        ti += blocks_tgt[k];
    }
    return out;
}

// Resize a captured dual onto the target operator's blocks and rescale by spec_norm_new / spec_norm_old:
// y lives in the spec_norm-normalized Ax-space and spec_norm depends on dt (and N), so this keeps the
// physical dual y / spec_norm unchanged.
inline Eigen::VectorXd ws_resize_dual(const OpWarmState &st, const std::vector<int> &blocks_tgt,
                                      double spec_norm_new) {
    Eigen::VectorXd y = ws_resize_blocks(st.y, st.blocks, blocks_tgt);
    if (st.spec_norm > 0.0) y *= spec_norm_new / st.spec_norm;
    return y;
}

// One maximal free/fixed run within a single axis-slice of a mask.
struct WsSeg {
    bool is_free;
    int start; // offset within the axis-slice
    int len;
};

// Split one axis-slice (mask[off .. off+n)) into maximal free(!=0)/fixed(==0) runs.
inline std::vector<WsSeg> ws_segments(const Eigen::VectorXd &mask, int off, int n) {
    std::vector<WsSeg> segs;
    int i = 0;
    while (i < n) {
        const bool is_free = mask(off + i) != 0.0;
        int j = i + 1;
        while (j < n && ((mask(off + j) != 0.0) == is_free)) j++;
        segs.push_back({is_free, i, j - i});
        i = j;
    }
    return segs;
}

// Resize a waveform run by run: free run k is resampled from source free run k (zeros if the source has
// fewer runs); fixed runs come from the target set_vals. Callers check free-run counts match first.
inline Eigen::VectorXd ws_resize_waveform(const Eigen::VectorXd &x_src, const Eigen::VectorXd &src_mask,
                                          const Eigen::VectorXd &tgt_mask, const Eigen::VectorXd &tgt_setvals,
                                          int Naxis) {
    const int n_src = static_cast<int>(src_mask.size()) / Naxis;
    const int n_tgt = static_cast<int>(tgt_mask.size()) / Naxis;
    Eigen::VectorXd out(tgt_mask.size());
    for (int ax = 0; ax < Naxis; ax++) {
        const int so = ax * n_src, to = ax * n_tgt;
        const std::vector<WsSeg> ssegs = ws_segments(src_mask, so, n_src);
        const std::vector<WsSeg> tsegs = ws_segments(tgt_mask, to, n_tgt);

        std::vector<Eigen::VectorXd> src_free; // source free blocks, in order
        for (const auto &s : ssegs) {
            if (s.is_free) src_free.push_back(x_src.segment(so + s.start, s.len));
        }

        size_t fi = 0;
        for (const auto &t : tsegs) {
            if (t.is_free) {
                out.segment(to + t.start, t.len) =
                    (fi < src_free.size()) ? ws_resample(src_free[fi], t.len) : Eigen::VectorXd::Zero(t.len);
                fi++;
            } else {
                out.segment(to + t.start, t.len) = tgt_setvals.segment(to + t.start, t.len);
            }
        }
    }
    return out;
}

// Count free runs on each axis (source and target must match before ws_resize_waveform).
inline std::vector<int> ws_free_run_counts(const Eigen::VectorXd &mask, int Naxis) {
    const int n = static_cast<int>(mask.size()) / Naxis;
    std::vector<int> counts(Naxis, 0);
    for (int ax = 0; ax < Naxis; ax++) {
        for (const auto &s : ws_segments(mask, ax * n, n)) {
            if (s.is_free) counts[ax]++;
        }
    }
    return counts;
}

} // namespace Gropt

#endif

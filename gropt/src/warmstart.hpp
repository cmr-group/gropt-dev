#ifndef GROPT_WARMSTART_H
#define GROPT_WARMSTART_H

// Warm-starting one solve from another (e.g. sweeping TE, or adding a constraint).
//
// What we carry, and why
//   * primal  X        -- the waveform. Resized across a grid change by interpolating only
//                         the FREE segments of the fixer mask; FIXED segments come from the
//                         new problem's set_vals (see ws_resize_waveform).
//   * dual    y (per op) -- the Lagrange multiplier. This is the high-value state: it
//                         accumulates over iterations and cannot be cheaply regenerated.
//   * weight rho, gamma  -- seeded; the reweighter (Xu et al.) re-adapts them.
// We deliberately do NOT carry the consensus z: it is regenerated as z = A*X after X is
// loaded (the cold-start path already does exactly this in WorkspaceSolver::prep).
//
// Operators are matched between solves by Operator::unique_name (assigned in prepare() as
// "<name>#<occurrence>"), so duplicate-named operators stay distinguishable and a rebuilt,
// resized operator set still maps correctly as long as it is built in the same order. An
// operator with no matching key (e.g. a newly added constraint) is left cold.

#include "Eigen/Dense"
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace Gropt {

// Snapshot of one operator's ADMM dual state. `blocks` records the flat partition of `y`
// at capture time (e.g. SAFE = n_terms*Naxis blocks of N) so the loader can resize `y`
// onto a new grid without needing the original operator object.
struct OpWarmState {
    std::string key;
    Eigen::VectorXd y;       // dual / Lagrange multiplier, in the operator's Ax-space
    double weight = 1.0;     // ADMM penalty rho_i (seed; reweighter re-adapts)
    double gamma = 1.5;      // relaxation gamma_i
    double spec_norm = 1.0;  // operator normalization at capture; renormalizes the dual on load (dt/N)
    std::vector<int> blocks; // partition of y at capture (sum == y.size())
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
// n_new == n returns a copy; degenerate sizes are handled gracefully.
inline Eigen::VectorXd ws_resample(const Eigen::VectorXd &v, int n_new) {
    const int n = static_cast<int>(v.size());
    if (n_new <= 0) return Eigen::VectorXd();
    if (n_new == n) return v;
    if (n <= 1) return Eigen::VectorXd::Constant(n_new, n == 1 ? v(0) : 0.0);
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

// Resize a vector partitioned into blocks_src to the partition blocks_tgt, interpolating
// each block independently so stacked blocks (e.g. SAFE's three term-blocks) never bleed
// into each other. If the block COUNTS differ the input is returned unchanged (the caller
// validates and warns -- a count change means the operator type changed, not just N).
inline Eigen::VectorXd ws_resize_blocks(const Eigen::VectorXd &v, const std::vector<int> &blocks_src,
                                        const std::vector<int> &blocks_tgt) {
    if (blocks_src.size() != blocks_tgt.size()) return v;
    int tgt_total = 0;
    for (int b : blocks_tgt) tgt_total += b;
    Eigen::VectorXd out(tgt_total);
    int si = 0, ti = 0;
    for (size_t k = 0; k < blocks_src.size(); k++) {
        out.segment(ti, blocks_tgt[k]) = ws_resample(v.segment(si, blocks_src[k]), blocks_tgt[k]);
        si += blocks_src[k];
        ti += blocks_tgt[k];
    }
    return out;
}

// Resize a captured dual onto a target operator's grid AND renormalize it for that operator's
// spec_norm. The snapshot's dual lives in the SOURCE operator's spec_norm-normalized space, and
// spec_norm depends on dt (and N for some operators), so preserving the *physical* dual across a
// grid change requires scaling by spec_norm_new / spec_norm_old after the block interpolation.
// blocks_tgt = target operator's Ax_block_lengths(); spec_norm_new = target operator's spec_norm.
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

// Segment-aware resize of a per-axis waveform. FREE runs are interpolated source-run k ->
// target-run k; FIXED runs are taken from the target problem's set_vals. Requires equal
// free-run counts per axis (the topology contract); missing source runs fill with zeros.
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

// Count free runs per axis (used to validate the topology contract before resizing X).
inline int ws_free_run_count(const Eigen::VectorXd &mask, int Naxis) {
    const int n = static_cast<int>(mask.size()) / Naxis;
    int count = 0;
    for (const auto &s : ws_segments(mask, 0, n)) {
        if (s.is_free) count++;
    }
    return count;
}

} // namespace Gropt

#endif

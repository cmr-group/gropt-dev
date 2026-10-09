// diff_init layouts at times that are exactly on the grid: the same TE must give the same layout however the
// double was computed (a search's T_readout + (m + 0.5) dt versus the same TE typed in).

#include "gropt_params.hpp"
#include "test_util.hpp"

using namespace Gropt;
using namespace gropt_test;

namespace {

enum class Kind { simple, deadtime, preencode };

GroptParams layout(Kind kind, double dt, double TE, double T_pre = 0.0) {
    GroptParams gp;
    if (kind == Kind::simple) gp.diff_init(dt, TE, 3e-3, 5e-3, 16e-3);
    else if (kind == Kind::deadtime) gp.diff_init_deadtime(dt, TE, 3e-3, 5e-3, 16e-3);
    else gp.diff_init_preencode(dt, TE, 3e-3, 5e-3, 16e-3, T_pre);
    return gp;
}

bool same_layout(const GroptParams &a, const GroptParams &b) {
    if (a.pdata.N != b.pdata.N) return false;
    const Eigen::ArrayXd sa = a.pdata.set_vals.array().isNaN().select(9.0, a.pdata.set_vals.array());
    const Eigen::ArrayXd sb = b.pdata.set_vals.array().isNaN().select(9.0, b.pdata.set_vals.array());
    return (sa == sb).all() && (a.pdata.inv_vec.array() == b.pdata.inv_vec.array()).all();
}

// Index of the last fixed sample of the 180 gap (the first free sample after it, minus one).
int gap_end(const GroptParams &gp) {
    int i = gp.pdata.N - 2;
    while (i > 0 && std::isnan(gp.pdata.set_vals(i))) i--;   // back over the free post-180 run
    return i;
}

} // namespace

int run_diff_init_tests() {
    std::printf("\n=== diff_init on-grid rounding ===\n");
    int fail = 0;

    // TE = 62.2 ms at 400 us: (TE/2 + T_180/2) / dt is exactly 84, so the 180 gap ends at sample 84.
    {
        const double dt = 400e-6, te_sum = 16e-3 + 115.5 * dt;   // 0.062200000000000005
        GroptParams a = layout(Kind::simple, dt, te_sum), b = layout(Kind::simple, dt, 0.0622);
        fail += report(same_layout(a, b), "TE 62.2 ms: computed and typed TE give one layout");
        fail += report(gap_end(b) == 84, "TE 62.2 ms: 180 gap ends at the exact sample 84", gap_end(b));
    }
    // T_pre = 10 ms at 10 us is exactly 1000 pre-encoding samples (10e-3 / 10e-6 = 999.9999999999999).
    {
        GroptParams gp;
        const int n_pre = gp.diff_init_preencode(10e-6, 80e-3, 3e-3, 5e-3, 16e-3, 10e-3);
        fail += report(n_pre == 1000, "preencode: T_pre 10 ms at 10 us is 1000 samples", n_pre);
    }
    // Every bin of every layout kind: the bin-mid TE as a search computes it and as printed to 1 ns agree.
    for (Kind kind : {Kind::simple, Kind::deadtime, Kind::preencode}) {
        int n_bad = 0, n = 0;
        for (double dt : {200e-6, 250e-6, 270e-6, 400e-6}) {
            for (int m = static_cast<int>(0.03 / dt); m < static_cast<int>(0.12 / dt); m++, n++) {
                const double te = 16e-3 + (m + 0.5) * dt;
                const double typed = std::round(te * 1e9) / 1e9;
                n_bad += !same_layout(layout(kind, dt, te, 2e-3), layout(kind, dt, typed, 2e-3));
            }
        }
        const char *name = kind == Kind::simple ? "diff_init" : kind == Kind::deadtime ? "diff_init_deadtime"
                                                                                        : "diff_init_preencode";
        fail += report(n_bad == 0, std::string(name) + ": bin-mid TE computed vs printed, " + std::to_string(n) +
                                       " bins, layouts differing",
                       n_bad);
    }
    return fail;
}

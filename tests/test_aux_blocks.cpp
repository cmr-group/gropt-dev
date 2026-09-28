// ProblemData auxiliary primal blocks: registered by name, so operators needing the same quantity (the
// PNS and cardiac Op_SAFE both want u >= |slew|) share one block.

#include "test_util.hpp"

using namespace Gropt;
using namespace gropt_test;

int run_aux_block_tests() {
    std::printf("\nProblemData aux blocks\n");
    int failures = 0;

    ProblemData p;
    p.N = 100;
    p.Naxis = 3;
    const int nw = 300;

    failures += report(p.n_wave() == nw, "n_wave = N * Naxis");
    failures += report(p.n_total() == nw, "no blocks declared: total primal is the waveform");
    failures += report(p.aux_offset("abs_slew") == -1, "undeclared block reports offset -1");

    const int off1 = p.add_aux("abs_slew", nw);
    failures += report(off1 == nw, "first block starts right after the waveform");
    failures += report(p.n_total() == 2 * nw, "total primal grew by the block size");
    failures += report(p.aux_offset("abs_slew") == nw, "declared block is findable by name");

    // A second operator wanting the same quantity shares it.
    const int off2 = p.add_aux("abs_slew", nw);
    failures += report(off2 == off1, "redeclaring the same name returns the same offset");
    failures += report(p.n_total() == 2 * nw, "redeclaring does not grow the primal");

    const int off3 = p.add_aux("other", 50);
    failures += report(off3 == 2 * nw, "a second distinct block packs after the first");
    failures += report(p.n_total() == 2 * nw + 50, "total primal accounts for both blocks");
    failures += report(p.aux_offset("abs_slew") == nw, "the first block keeps its offset");

    // Disagreeing about a shared block's size is a programming error, not a resize.
    failures += report(throws_invalid([&] { p.add_aux("abs_slew", nw + 1); }),
                       "redeclaring a name with a different size throws");
    failures += report(p.n_total() == 2 * nw + 50, "the failed redeclaration left the primal unchanged");
    failures += report(throws_invalid([&] { p.add_aux("empty", 0); }), "a zero-sized block throws");

    // prepare() rebuilds the registry every call, so repeated prepares must not accumulate blocks.
    p.clear_aux();
    failures += report(p.n_total() == nw && p.aux.empty(), "clear_aux returns the primal to the waveform");
    failures += report(p.add_aux("abs_slew", nw) == nw, "re-declaring after a clear reuses the offset");

    return failures;
}

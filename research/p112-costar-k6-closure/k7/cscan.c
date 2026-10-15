#include <stdint.h>
/* For each clique q (bits split into lo/hi 64-bit words) and groups g (R,C masks):
   ar = |q&R_g|, ac = |q&C_g|; tot = sum_g min(cap_g, ar, ac); ct = sum_{g: ar>0 && ac>0} cap_g.
   viol: argmin over {tot < sz} of (tot - sz), first index on ties (-1 if none).
   tight: indices with ct == sz (in index order) with gamma bitmask; returns count (capped at maxt). */
int scan(int nq, const uint64_t *qlo, const uint64_t *qhi, const int8_t *qsz,
         int G, const uint64_t *Rlo, const uint64_t *Rhi, const uint64_t *Clo, const uint64_t *Chi,
         const int *caps, int *viol, int *tidx, uint32_t *tgam, int maxt)
{
    int best = 1000, bi = -1, nt = 0;
    for (int i = 0; i < nq; i++) {
        uint64_t lo = qlo[i], hi = qhi[i];
        int tot = 0, ct = 0; uint32_t gm = 0;
        for (int g = 0; g < G; g++) {
            int ar = __builtin_popcountll(lo & Rlo[g]) + __builtin_popcountll(hi & Rhi[g]);
            if (!ar) continue;
            int ac = __builtin_popcountll(lo & Clo[g]) + __builtin_popcountll(hi & Chi[g]);
            if (!ac) continue;
            int m = ar < ac ? ar : ac; if (caps[g] < m) m = caps[g];
            tot += m; ct += caps[g]; gm |= 1u << g;
        }
        int s = qsz[i];
        if (tot < s && tot - s < best) { best = tot - s; bi = i; }
        if (ct == s) { if (nt < maxt) { tidx[nt] = i; tgam[nt] = gm; } nt++; }
    }
    *viol = bi;
    return nt;
}

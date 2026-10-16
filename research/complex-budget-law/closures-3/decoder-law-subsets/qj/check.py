"""Independent checker (numpy; does not use sub.c / cand.c / bb*.c).
Subspaces are parametrised through the *form* (X = {u : B(u, g) = 0 for g in G}), not RREF of the
standard dual; duplicates are allowed (only maxima/minima are taken).
Checks: (1) thresholds used by the proof for all 255 support classes; (2) the low-root bound
max sum(delta) over simultaneously low root lines, for p <= PMAX_CHECK, by a vectorised exact B&B."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys, numpy as np, e8m
NL = 120
PMAX_CHECK = min(36, int(sys.argv[1]) if len(sys.argv) > 1 else 36)
def tomask(m): return np.array([m & ((1 << 64) - 1), m >> 64], dtype=np.uint64)
# hyperplanes via the form: lines r with B(u_r, g) = 0
Hm = np.zeros((512, 2), dtype=np.uint64)
for g in range(1, 512):
    m = sum(1 << r for r in range(NL) if e8m.form(e8m.LAB[r], g) == 0); Hm[g] = tomask(m)
G1, G2 = np.triu_indices(511, 1); G1 = G1 + 1; G2 = G2 + 1          # functionals g1 < g2 (distinct => independent)
pairs = Hm[G1] & Hm[G2]                                              # all codim-2 subspaces (each 3 times)
Hm = Hm[1:]
def pc(A): return (np.bitwise_count(A[..., 0]) + np.bitwise_count(A[..., 1])).astype(np.int64)
ok = True
# (1) thresholds
res = {}
for nodes, S in e8m.CLASSES:
    s = tomask(S); n = bin(S).count('1')
    t8 = n - pc(Hm & s).max(); t7 = n - pc(pairs & s).max()
    res.setdefault(n, []).append((t8, t7))
root = res[57]; d7 = res[78]
other_t7 = min(t7 for n, v in res.items() if n not in (57, 78) for t8, t7 in v)
print("root (to<=8, to<=7):", set(root), " D7:", set(d7), " other classes min to<=7:", other_t7, flush=True)
ok &= set(root) == {(1, 25)} and set(d7) == {(14, 42)} and other_t7 == 51
# root span<=6 and <=5 thresholds (codim 3, 4; with repeats) for line 0
N0 = e8m.supp(e8m.LINES[0]); s0 = tomask(N0)
t6 = 99; t5 = 99
for g3 in range(1, 512):
    keep = (G1 != g3) & (G2 != g3) & ((G1 ^ G2) != g3)                 # g3 independent of g1, g2 => codim 3
    tri = pairs[keep] & Hm[g3 - 1]; t6 = min(t6, 57 - pc(tri & s0).max())
print("root to<=6:", t6, flush=True); ok &= (t6 == 37)   # span<=5 (45) is checked by sub.c only
# (2) witnesses (span<=7 and span<=6) for every root line, size <= PMAX_CHECK
W, Dl, Own = [], [], []
for a in range(NL):
    Na = tomask(e8m.supp(e8m.LINES[a]))
    w7 = Na & ~pairs; w7 = np.unique(w7[pc(w7) <= PMAX_CHECK], axis=0)
    w6 = np.zeros((0, 2), np.uint64)   # span<=6 needs >= 37 deletions (checked above): irrelevant for p <= 36
    W += [w7, w6]; Dl += [np.ones(len(w7), int), 2 * np.ones(len(w6), int)]; Own += [np.full(len(w7), a), np.full(len(w6), a)]
W = np.concatenate(W); Dl = np.concatenate(Dl); Own = np.concatenate(Own)
print("witnesses:", len(W), " root0:", int((Own == 0).sum()))
def maxdelta(P):
    sel = pc(W) <= P; Wp, Dp, Op = W[sel], Dl[sel], Own[sel]
    best = 0; stack = [(int(i), Wp[i], int(Dp[i])) for i in np.nonzero(Op == 0)[0]]
    nodes = 0
    while stack:
        i, U, d = stack.pop(); nodes += 1; best = max(best, d)
        cand = np.nonzero((Op > Op[i]) & (pc(Wp | U) <= P))[0]
        for j in cand: stack.append((int(j), Wp[j] | U, d + int(Dp[j])))
    return best, nodes
expect = {25: 1, 26: 2, 30: 2, 35: 2, 36: 6, 37: 6, 38: 6, 39: 6, 40: 6, 41: 6}
for P in range(25, PMAX_CHECK + 1):
    b, nodes = maxdelta(P); lim = ((P - 8) // 2) // 2
    good = b <= lim and (P not in expect or expect[P] == b)
    print("p=%d max sum delta=%d limit=%d nodes=%d %s" % (P, b, lim, nodes, "PASS" if good else "FAIL"), flush=True)
    ok &= good
print("CHECK", "PASS" if ok else "FAIL")
print("PASS_EXPQJ_CHECK" if ok else "FAIL_EXPQJ_CHECK")

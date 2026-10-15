"""Independent k=7 orbit completeness: every rep is the lex-least sorted image of itself under ALL of EL
(brute force, no shortcut) => reps pairwise inequivalent; own stabiliser counts; sum |EL|/|Stab| == C(120,7).
Uses only numpy, group_elems.npy and reps7.npz (no model/engine/stage code).
usage: orb_indep.py DIR [JOBS [START STOP]]   DIR holds group_elems.npy and reps7.npz (e.g. "."); JOBS default 3.
With START STOP only reps[START:STOP] (covering reps first, then non-covering) are checked and the sum test
is skipped (smoke test)."""
import numpy as np, sys, math
from multiprocessing import Pool
AP = sys.argv[1]; JOBS = int(sys.argv[2]) if len(sys.argv) > 2 else 3
SL = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else None
EL = np.load(AP + '/group_elems.npy').astype(np.int64); N = len(EL)
z = np.load(AP + '/reps7.npz')
dec = lambda c: [(c >> (7 * (6 - i))) & 127 for i in range(7)]
reps = [list(map(int, a)) for a in z['cov_A']] + [dec(int(c)) for c in z['non_code']]
if SL: reps = reps[SL[0]:SL[1]]
W = 128 ** np.arange(6, -1, -1, dtype=np.int64)
def work(chunk):
    out = []
    for A in chunk:
        img = np.sort(EL[:, A], axis=1) @ W
        a = np.array(sorted(A)) @ W
        out.append((int(img.min() == a and A == sorted(A)), int((img == a).sum())))
    return out
if __name__ == '__main__':
    ch = [reps[i:i + 2000] for i in range(0, len(reps), 2000)]
    with Pool(JOBS) as p: res = [r for c in p.imap(work, ch) for r in c]
    assert len(set(map(tuple, reps))) == len(reps)
    print('reps', len(reps), 'all lex-least', all(m for m, _ in res))
    if SL: print('slice', SL, 'stab sizes divide |EL|', all(N % s == 0 for _, s in res)); sys.exit()
    tot = sum(N // s for _, s in res); assert all(N % s == 0 for _, s in res)
    print('sum', tot, 'C(120,7)', math.comb(120, 7), 'MATCH' if tot == math.comb(120, 7) else 'MISMATCH')

"""Tally stage-2 results per composition shape; write survivors to cases_k7_all.json."""
import glob, pickle, json, numpy as np, itertools
from collections import Counter, defaultdict
from model import GM
COMPS = [c for c in itertools.product(range(1, 10), repeat=7) if sum(c) == 9]
NAMES = ["FILTER_L1", "L0", "L1", "L1Z", "L3odd", "L2rank", "L3GZodd", "SURVIVE", "ERROR"]
z = np.load("reps7.npz"); nrep = len(z["cov_code"])
ds = sorted(glob.glob("s2/d*.npz")); assert len(ds) == len(glob.glob("s2/d*.done"))
R = np.concatenate([np.load(f)["rows"] for f in ds])
print("chunks", len(ds), "rows", len(R), "expected", nrep * 28)
reps_seen = np.unique(R[:, 0]); print("reps covered", len(reps_seen), "of", nrep)
cnt = defaultdict(Counter)
for ci in range(28):
    shape = "(3,1^6)" if 3 in COMPS[ci] else "(2,2,1^5)"
    sub = R[R[:, 1] == ci]
    for k, n in zip(*np.unique(sub[:, 2], return_counts=True)): cnt[shape][NAMES[k]] += int(n)
for s in cnt: print(s, dict(cnt[s]), "total", sum(cnt[s].values()))
closed = R[(R[:, 2] < 7)]
print("certificates", len(closed), "replay ok", int(closed[:, 4].sum()), "FAIL", int((closed[:, 4] == 0).sum()))
surv = []
for f in sorted(glob.glob("s2/v*.pkl")):
    for A, c, fin in pickle.load(open(f, "rb")): surv.append((list(A), list(c)))
t = Counter()
for A, c in surv:
    big = [A[i] for i in range(7) if c[i] > 1]
    t[(tuple(sorted(c, reverse=True)), "A2pair" if len(big) == 2 and GM[big[0]][big[1]] != 0 else "orthpair" if len(big) == 2 else "-")] += 1
print("survivors", len(surv), dict(t))
json.dump(surv, open("cases_k7_all.json", "w"))

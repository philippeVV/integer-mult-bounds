"""Merge stage-1 chunks: one representative per covering orbit (first in chunk order), non-covering codes."""
import glob, numpy as np, sys
import os
fs = sorted(glob.glob("s1/c*.npz")); PART = os.environ.get("PARTIAL")
if not PART: assert len(fs) == 323 and len(glob.glob("s1/c*.done")) == 323
else: fs = [f for f in fs if os.path.exists(f[:-4] + ".done")]
cov = {}; non = set(); cnt = np.zeros(4, dtype=np.int64)
for f in fs:
    z = np.load(f); cnt += z["cnt"]
    for k, A, m in zip(z["cov_code"].tolist(), z["cov_A"].tolist(), z["cov_mask"].tolist()):
        if k in cov:
            assert bin(cov[k][1]).count("1") == bin(m).count("1"), k
        else:
            cov[k] = (A, m)
    non.update(z["non_code"].tolist())
assert not (non & set(cov)), "covering and non-covering codes overlap"
ks = sorted(cov)
A = np.array([cov[k][0] for k in ks], dtype=np.uint8); M = np.array([cov[k][1] for k in ks], dtype=np.int32)
np.savez_compressed(os.environ.get("REPS", "reps7.npz"), cov_code=np.array(ks, dtype=np.int64), cov_A=A, cov_mask=M,
                    non_code=np.array(sorted(non), dtype=np.int64))
pc = np.array([bin(m).count("1") for m in M.tolist()])
print("extension counts noncov/cov/L1kills/survpairs(raw):", cnt.tolist())
print("covering orbits", len(ks), "non-covering orbits", len(non), "total 7-orbits", len(ks) + len(non))
print("covering orbits with >=1 filter survivor", int((pc > 0).sum()), "surviving (orbit, comp) pairs", int(pc.sum()))

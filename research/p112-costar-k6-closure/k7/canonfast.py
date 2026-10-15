"""Faster canonical form for subsets of lines under H.
Canonical image = lexicographically least sorted image.  Its least element is 0 if the
set meets F0 (F0 is one H-orbit containing 0) and 2 otherwise (live lines form one
H-orbit containing 2).  Hence only group elements sending some member of the
relevant orbit to that minimum need be scanned."""
import numpy as np
from model import F0

EL = np.load("group_elems.npy").astype(np.int16)   # (172032,120)
Wt = np.array([128 ** i for i in range(9)][::-1], dtype=np.int64)
F0S = set(F0)
TO0 = {a: EL[EL[:, a] == 0] for a in F0}
TO2 = {a: EL[EL[:, a] == 2] for a in range(120) if a not in F0S}
assert all(len(v) == 21504 for v in TO0.values()) and all(len(v) == 1536 for v in TO2.values())

def canon(A):
    A = list(A)
    k = len(A)
    f = [a for a in A if a in F0S]
    if f:
        cand = np.concatenate([TO0[a] for a in f])
    else:
        cand = np.concatenate([TO2[a] for a in A])
    img = cand[:, A].astype(np.int64)
    img.sort(axis=1)
    code = img @ Wt[-k:]
    i = int(np.argmin(code))
    return tuple(int(x) for x in img[i])

def stab_size(A):
    A = sorted(A)
    k = len(A)
    img = EL[:, A].astype(np.int64)
    img.sort(axis=1)
    code = img @ Wt[-k:]
    c0 = np.array(A, dtype=np.int64) @ Wt[-k:]
    return int((code == c0).sum())

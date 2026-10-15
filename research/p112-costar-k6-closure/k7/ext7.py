"""k=7 initial filter, used by stage1.py and stage2.py.
For a 6-orbit rep B (../engine/orbk.pkl[6]) and x notin B, A = B+{x}:
  not covering -> L0 at full supports (some live line in no co-star): infeasible for every composition
  covering     -> initial clique filter (L1 at full co-star supports) per composition of 9 into 7 parts.
filt(A) and filt_fast(B) return the surviving compositions and, for each killed one, the index of a
violating clique in cliques.pkl.  (The prototype driver that used to live here is removed; filt,
filt_fast and init are unchanged from the run.)"""
import pickle, sys, time, os, random, itertools
import numpy as np
from multiprocessing import Pool
from model import COSTAR, FULL, NL, bits
COMPS = [c for c in itertools.product(range(1, 10), repeat=7) if sum(c) == 9]
_MQ = _SZ = _CL = None
def init():
    global _MQ, _SZ, _CL
    _CL = pickle.load(open("cliques.pkl", "rb"))
    M = np.zeros((len(_CL), NL), dtype=np.int16)
    for i, q in enumerate(_CL):
        for b in bits(q): M[i, b] = 1
    V = np.zeros((NL, 120), dtype=np.int16)
    for a in range(120):
        for b in bits(COSTAR[a]): V[b, a] = 1
    _MQ = (M @ V).astype(np.int8)
    _SZ = M.sum(axis=1).astype(np.int16)
def filt(A):
    MA = _MQ[:, list(A)]
    T1 = (MA > 0).sum(axis=1, dtype=np.int16)
    dq = np.nonzero(T1 < _SZ)[0]
    MAd = MA[dq].astype(np.int16); SZd = _SZ[dq]
    surv = []; dead = {}
    for c in COMPS:
        tot = np.minimum(MAd, np.array(c, dtype=np.int16)).sum(axis=1)
        bad = np.nonzero(tot < SZd)[0]
        if len(bad): dead[c] = int(dq[bad[0]])
        else: surv.append(c)
    return surv, dead
_P = _B2 = _B3 = None
C2 = np.array([[1 if c[g] >= 2 else 0 for c in COMPS] for g in range(7)], dtype=np.int16)
C3 = np.array([[1 if c[g] >= 3 else 0 for c in COMPS] for g in range(7)], dtype=np.int16)
def filt_fast(B):
    """yields (x, A, surv, dead) for all x notin B, A covering"""
    global _P, _B2, _B3
    if _P is None:
        _P = (_MQ > 0).astype(np.int8); _B2 = (_MQ >= 2).astype(np.int8); _B3 = (_MQ >= 3).astype(np.int8)
    TB = _P[:, list(B)].sum(axis=1, dtype=np.int16)
    D0 = _SZ - TB                      # deficit before adding x
    UB = 0
    for a in B: UB |= COSTAR[a]
    cand = np.nonzero(D0 > 0)[0]
    Pc = _P[cand]                      # (nc,120)
    for x in range(120):
        if x in B: continue
        if (UB | COSTAR[x]) != FULL:
            yield x, None, None, None; continue
        A = tuple(sorted(B + (x,)))
        dl = D0[cand] - Pc[:, x]
        sel = np.nonzero(dl > 0)[0]
        if len(sel) == 0:
            yield x, A, list(COMPS), {}; continue
        qi = cand[sel]; dd = dl[sel]
        b2 = _B2[np.ix_(qi, list(A))].astype(np.int16); b3 = _B3[np.ix_(qi, list(A))].astype(np.int16)
        gain = b2 @ C2 + b3 @ C3 - dd[:, None]      # (nd,28)
        ok = (gain >= 0).all(axis=0)
        surv = [COMPS[i] for i in np.nonzero(ok)[0]]
        dead = {}
        for i in np.nonzero(~ok)[0]:
            dead[COMPS[i]] = int(qi[int(np.argmax(gain[:, i] < 0))])
        yield x, A, surv, dead

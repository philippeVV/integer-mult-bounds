"""Sub-configurations L of the 120 E8 lines: regular subsystems (r.v = 0 mod p), flats (orthogonal to roots)."""
from e8 import *
import itertools
def lattice_vectors(maxnorm2=12):
    # E8 vectors in doubled coords: all-even-integer coords (x2) with even sum/... use sums of <=3 roots
    R = [r for r in ROOT] + [tuple(-x for x in r) for r in ROOT]
    vs = set(R)
    for a in R:
        for b in R:
            vs.add(tuple(x + y for x, y in zip(a, b)))
    return [v for v in vs if any(v)]
def mod_sub(v, p):
    """lines r with (r.v) = 0 mod p, true inner product r.v = dot/4"""
    return sum(1 << a for a in range(NL) if (sum(x * y for x, y in zip(ROOT[a], v)) // 4) % p == 0)
def flat(roots_idx):
    """lines orthogonal to all given lines"""
    m = FULL
    for i in roots_idx: m &= ORTH[i]
    return m
def components(L):
    """connected components of non-orthogonality graph on L (= irreducible root subsystems)"""
    comps = []; rest = L
    while rest:
        a = (rest & -rest).bit_length() - 1; comp = 1 << a; front = comp
        while front:
            b = (front & -front).bit_length() - 1; front ^= 1 << b
            nb = (L & ~ORTH[b] & ~comp) & ~(1 << b)
            comp |= nb; front |= nb
        comps.append(comp); rest &= ~comp
    return comps
def typ(L):
    # identify by (lines, rank) per component
    out = []
    for c in components(L):
        k = bin(c).count('1')
        rk = rank_q([ROOT[a] for a in bits(c)])
        out.append((k, rk))
    return sorted(out, reverse=True)
def rank_q(vs):
    from fractions import Fraction as Fr
    A = [[Fr(x) for x in v] for v in vs]; r = 0
    for c in range(8):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c] / A[r][c]; A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        r += 1
    return r
NAMES = {(1,1):'A1',(3,2):'A2',(6,3):'A3',(10,4):'A4',(12,4):'D4',(15,5):'A5',(20,5):'D5',(21,6):'A6',(30,6):'D6',(36,6):'E6',
         (28,7):'A7',(42,7):'D7',(63,7):'E7',(36,8):'A8',(56,8):'D8',(120,8):'E8'}
def name(L): return '+'.join(NAMES.get(t, str(t)) for t in typ(L))
if __name__ == '__main__':
    found = {}
    V = lattice_vectors()
    for v in V:
        for p in (2, 3, 4, 5, 6):
            L = mod_sub(v, p)
            if L and L != FULL:
                nm = name(L)
                if rank_q([ROOT[a] for a in bits(L)]) == 8 and nm not in found: found[nm] = L
    for nm, L in sorted(found.items(), key=lambda x: -bin(x[1]).count('1')): print(nm, bin(L).count('1'))

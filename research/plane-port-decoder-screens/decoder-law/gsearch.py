"""Anneal over lattice-gauge decoders N = G/2 restricted to L (basis M of 8 rays, dual basis), maximise
gain(L, M) = |L| - sum_k rho(B[R_k cap L, A_k cap L]).  Exact (G/2 restricted to any L is an exact masked decoder)."""
import sys, random, math
from e8 import *
from fractions import Fraction as Fr
def nullvec(rows):
    """integer vector orthogonal to the given 7 integer vectors (or None if rank < 7)."""
    n = 8; A = [[Fr(x) for x in r] for r in rows]; piv = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]; A[r] = [x / A[r][c] for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    if r != 7: return None
    free = [c for c in range(n) if c not in piv][0]
    v = [Fr(0)] * n; v[free] = Fr(1)
    for i, c in enumerate(piv): v[c] = -A[i][free]
    from math import lcm
    d = 1
    for x in v: d = lcm(d, x.denominator)
    return tuple(int(x * d) for x in v)
def rand_ray(rnd):
    while True:
        v = nullvec([ROOT[rnd.randrange(NL)] for _ in range(7)])
        if v: return v
def indep(M):
    try: inv_dual(M); return True
    except StopIteration: return False
def gain(M, L):
    return bin(L).count('1') - gauge_cost(M, L)
if __name__ == '__main__':
    seed, iters = int(sys.argv[1]), int(sys.argv[2]); start = sys.argv[3] if len(sys.argv) > 3 else 'frame'
    rnd = random.Random(seed)
    if start == 'frame': M = [ROOT[a] for a in frames()[rnd.randrange(2025)]]
    else:
        while True:
            M = [rand_ray(rnd) for _ in range(8)]
            if indep(M): break
    L = FULL; cur = gain(M, L); best = cur; T = 1.5
    for it in range(iters):
        if rnd.random() < 0.3:
            k = rnd.randrange(8); old = M[k]; M[k] = rand_ray(rnd) if rnd.random() < 0.7 else ROOT[rnd.randrange(NL)]
            if not indep(M): M[k] = old; continue
            new = gain(M, L)
            if new >= cur or rnd.random() < math.exp((new - cur) / T): cur = new
            else: M[k] = old
        else:
            a = rnd.randrange(NL); L2 = L ^ (1 << a); new = gain(M, L2)
            if new >= cur or rnd.random() < math.exp((new - cur) / T): cur, L = new, L2
        T = max(0.1, T * 0.9995)
        if cur > best:
            best = cur; print('it', it, 'gain', best, '|L|', bin(L).count('1'), 'cost', bin(L).count('1') - best, flush=True)
    print('final best gain', best)

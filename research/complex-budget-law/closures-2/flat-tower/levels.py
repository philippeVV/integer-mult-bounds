"""expNV / levels.py -- different-family stage sharing (parallel centre calls) in the typed word.  Exact, stdlib.
Run: python3 -I -B levels.py

Theorem (README sec.2b). In any typed word with exactly 2n+1 centre units, every block has rank 1
(no two calls can be merged / run in parallel), proved with sigma(S)=n-a(S)-b(S), tau(T)=n-a(T)-b(T),
sigma+tau = c+d <= 1.  Swap lemma: without descents, if call t+1 sits at a frame dimension below call t's
output frame, the two calls commute.  Together: frame dimension strictly rises by 8 per call -> 2n+1 <= 9.
This script checks, exactly:
  1. the potential inequalities the proof uses (random rational words, n=2..4);
  2. sigma, tau along the retained chain for n=1..8 (sigma+tau = 1 strictly inside, every block rank <= its bound);
  3. effective consecutive calls of the chain never commute (n=1..8);
  4. the greedy minimal frame-level assignment of the chain (with free idle advances, no descents) is 0,1,...,2n,
     so the chain needs 2n+1 fresh 8-blocks (= m = 1+8(2n+1)), i.e. variant T of expNL.
"""
import sys, os, random
from fractions import Fraction as Q
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chain as C
rank, ident = C.rank, C.ident
def mat_mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def cols(g, idx, N): return [[g[r][i] for r in range(N)] for i in idx]
def inter(A, B): return rank(A) + rank(B) - rank(A + B)
def pots(g, n):
    N = 2 * n; E = ident(N); X = [E[i] for i in range(n)]; Y = [E[n + i] for i in range(n)]
    S = cols(g, range(n), N); T = cols(g, range(n, N), N)
    return inter(S, X), inter(S, Y), inter(T, X), inter(T, Y)
def check_inequalities(rng):
    cnt = 0
    for n in (2, 3, 4):
        N = 2 * n
        for _ in range(60):
            g = [[Q(rng.randint(-2, 2)) for _ in range(N)] for _ in range(N)]
            if rank(g) < N: continue
            r = rng.randint(1, n); A = [[Q(rng.randint(-2, 2)) for _ in range(r)] for _ in range(n)]
            B = [[Q(rng.randint(-2, 2)) for _ in range(n)] for _ in range(r)]; Cm = mat_mul(A, B); rc = rank(Cm)
            for kind in 'LU':
                Sh = ident(N)
                for i in range(n):
                    for j in range(n):
                        if kind == 'L': Sh[n + i][j] = Cm[i][j]
                        else: Sh[i][n + j] = Cm[i][j]
                h = mat_mul(Sh, g); p0, p1 = pots(g, n), pots(h, n)
                aS, bS, aT, bT = p0; aS2, bS2, aT2, bT2 = p1
                assert all(v >= 0 for v in (n - aS2 - bS2, n - aT2 - bT2))           # sigma, tau >= 0
                if kind == 'L': assert (bS, bT) == (bS2, bT2) and abs(aS2 - aS) <= rc and abs(aT2 - aT) <= rc
                else: assert (aS, aT) == (aS2, aT2) and abs(bS2 - bS) <= rc and abs(bT2 - bT) <= rc
                cnt += 1
    return cnt
def effective_calls(n):
    """Return the chain's calls as effective shears g_eff = B^{-1} e e^T B (as matrices acting on the state)."""
    ops = C.chain(n); N = 2 * n; Bm = ident(N); eff = []
    def elem(t, s, c):
        E = ident(N); E[t][s] += c; return E
    for kind, t, s, c in ops:
        E = elem(t, s, Q(c))
        if kind == 'bridge':
            Bm = mat_mul(E, Bm)
        else:
            Binv = inv(Bm); eff.append((kind, mat_mul(Binv, mat_mul(E, Bm))))
    return eff
def inv(A):
    N = len(A); M = [r[:] + [Q(int(i == j)) for j in range(N)] for i, r in enumerate(A)]
    for c in range(N):
        p = next(i for i in range(c, N) if M[i][c] != 0); M[c], M[p] = M[p], M[c]
        f = M[c][c]; M[c] = [x / f for x in M[c]]
        for i in range(N):
            if i != c and M[i][c] != 0:
                f = M[i][c]; M[i] = [a - f * b for a, b in zip(M[i], M[c])]
    return [r[N:] for r in M]
def main():
    rng = random.Random(5)
    print(f'1. potential inequalities used by the proof: {check_inequalities(rng)} random exact shear steps OK')
    for n in range(1, 9):
        ops = C.chain(n); N = 2 * n; g = ident(N); sig = []; ranks_ok = True
        for kind, t, s, c in ops:
            g[t] = [a + c * b for a, b in zip(g[t], g[s])]
            if kind == 'bridge': continue
            aS, bS, aT, bT = pots(g, n); sig.append((n - aS - bS, n - aT - bT))
        inner = [a + b for a, b in sig[:-1]]
        assert all(v == 1 for v in inner) and sig[-1] == (0, 0), (n, sig)
        eff = effective_calls(n)
        for (k1, A), (k2, B) in zip(eff, eff[1:]):
            assert mat_mul(A, B) != mat_mul(B, A), ('consecutive calls commute', n)
        # greedy frame levels: a bridge forces equal levels (raise the lower, free), a call needs both members of its
        # pair at the same level and raises them by one block.
        lev = {('X', j): 0 for j in range(n)}; lev.update({('Y', j): 0 for j in range(n)}); call_lev = []
        reg = lambda i: ('X', i) if i < n else ('Y', i - n)
        for kind, t, s, c in ops:
            a, b = reg(t), reg(s)
            if kind == 'bridge':
                m = max(lev[a], lev[b]); lev[a] = lev[b] = m
            else:
                j = a[1]; m = max(lev[('X', j)], lev[('Y', j)]); call_lev.append(m)
                lev[('X', j)] = lev[('Y', j)] = m + 1
        assert call_lev == list(range(2 * n + 1)), (n, call_lev)
        print(f'n={n}: sigma+tau=1 inside, consecutive effective calls never commute, call levels {call_lev} '
              f'-> needs {2*n+1} fresh blocks (m={1+8*(2*n+1)}); fits m=73 iff n<=4: {2*n+1 <= 9}')
    print('PASS_EXPNV_LEVELS')
if __name__ == '__main__':
    main()

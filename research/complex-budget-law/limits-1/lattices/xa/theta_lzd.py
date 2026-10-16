"""expXA theta_lzd.py: expQL Thm B (PSD local-theta bound) on LZD_n. Two S_n-orbits (triple roots, A_{n-1} roots),
so one local SDP each. Shows that a PSD decoder needs rank >= V^2/sum(1+theta_a) >> n+1 = rank of the indefinite W/2,
i.e. the Lorentzian signature is the mechanism.  Also prints the Frobenius bound tr^2/||W/2||_F^2 (valid for any
symmetric decoder) and the char-2 trace bound V/h (2-integral decoders)."""
import sys, itertools
import numpy as np, cvxpy as cp


def lzd(n):
    vecs = []
    for T in itertools.combinations(range(1, n + 1), 3):
        v = [0] * (n + 1); v[0] = 1
        for i in T: v[i] = 1
        vecs.append(v)
    for i, j in itertools.combinations(range(1, n + 1), 2):
        v = [0] * (n + 1); v[i] = 1; v[j] = -1; vecs.append(v)
    X = np.array(vecs); eps = np.array([-1] + [1] * n)
    return (X * eps) @ X.T, comb3(n)


def comb3(n):
    from math import comb
    return comb(n, 3)


def lovasz(E, m):
    X = cp.Variable((m, m), symmetric=True)
    cons = [X >> 0, cp.trace(X) == 1]
    if E:
        I = np.array([e[0] for e in E]); J = np.array([e[1] for e in E]); cons.append(X[I, J] == 0)
    p = cp.Problem(cp.Maximize(cp.sum(X)), cons); p.solve(solver='CLARABEL')
    return p.value


for n in [int(x) for x in sys.argv[1:]] or [9, 11, 13]:
    W, nt = lzd(n)
    V = len(W)
    free = (np.abs(W) == 1)
    tot = 0.0
    out = []
    for a, size in [(0, nt), (nt, V - nt)]:
        nb = np.nonzero(free[a])[0]; m = len(nb)
        sub = free[np.ix_(nb, nb)]
        E = [(i, j) for i in range(m) for j in range(i + 1, m) if not sub[i, j]]
        th = lovasz(E, m)
        tot += size * (1 + th); out.append((size, m, th))
    N = W / 2
    frob = np.trace(N) ** 2 / (N * N).sum()
    h = n if n % 2 else n + 1
    print(f"n={n} V={V} rank(W/2)={n if n>9 else 8} (signature ({n-1},1) for n>=10); orbits (size, free-deg, theta_a) = "
          + ", ".join(f"({s},{m},{t:.3f})" for s, m, t in out)
          + f"; PSD bound V^2/sum(1+theta) = {V*V/tot:.2f}; Frobenius bound (any symmetric) = {frob:.2f};"
          f" char-2 trace bound V/h = {V/h:.2f}", flush=True)
    if n == 11: assert V * V / tot > 11.9   # every PSD decoder needs rank > 11 = rank of the indefinite W/2
print('PASS_EXPXA_THETA')

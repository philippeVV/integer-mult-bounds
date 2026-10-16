#!/usr/bin/env python3
"""expNQ: numerical real min-rank of the masked identity for orthogonality masks, and the gamma sandwich.

Mask: N[a,a] = 1, N[a,b] = 0 when ports a,b are NOT orthogonal (required), free when orthogonal.
Families:
  e8        : 120 E8 root lines, labels (r mod 2, 1) in F2^9; free iff <r,r'> odd   (theta = 1/9)
  pg n q    : F_q-lines (q = 2^u) in F_q^n = F2^{nu}, form sum x_i y_i (trace form after self-dual basis);
              ports = non-isotropic points of PG(n-1,q); free iff B(p,r) = 0          (theta = 1/n)
Sandwich (field-line labels, PR in {0,u}):  mr/V <= gamma_min <= mr/(theta V)   (README sec. 2)
F2 floor (Theorem A form):  rank_F2 of every F2 masked identity >= theta*V.
Method: alternating least squares on N = X Y^T restricted to required entries, many restarts, residual test.
usage: minrank.py e8 | minrank.py pg n q   [rmin rmax]
"""
import sys, itertools, numpy as np

def gf_tables(u):
    polys = {1: 0b11, 2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101}
    q = 1 << u; P = polys[u]
    def mul(a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            b >>= 1; a <<= 1
            if a & q: a ^= P
        return r
    M = np.array([[mul(a, b) for b in range(q)] for a in range(q)], dtype=np.int64)
    inv = [0] * q
    for a in range(1, q):
        for b in range(1, q):
            if M[a, b] == 1: inv[a] = b
    return M, inv

def pg_family(n, q):
    u = q.bit_length() - 1
    M, inv = gf_tables(u)
    pts = []
    for v in itertools.product(range(q), repeat=n):
        if not any(v): continue
        lead = next(x for x in v if x)
        if lead != 1: continue          # normalised representative
        pts.append(v)
    def B(x, y):
        s = 0
        for a, b in zip(x, y): s ^= int(M[a, b])
        return s
    pts = [p for p in pts if B(p, p) != 0]
    V = len(pts)
    free = np.zeros((V, V), bool)
    for i in range(V):
        for j in range(i + 1, V):
            if B(pts[i], pts[j]) == 0: free[i, j] = free[j, i] = True
    return free, 1.0 / n, n * u, f"PG({n-1},{q}) F{q}-lines, u={u}, h={n*u}"

def e8_family():
    roots = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = [0] * 8; v[i] = 2 * si; v[j] = 2 * sj; roots.append(tuple(v))   # doubled coords
    for s in itertools.product((1, -1), repeat=8):
        if s.count(-1) % 2 == 0: roots.append(tuple(s))
    lines = [r for r in roots if r > tuple(-x for x in r)]
    V = len(lines)
    ip = lambda a, b: sum(x * y for x, y in zip(a, b)) // 4
    free = np.zeros((V, V), bool)
    for i in range(V):
        for j in range(i + 1, V):
            if ip(lines[i], lines[j]) % 2: free[i, j] = free[j, i] = True
    return free, 1.0 / 9, 9, "E8 P120 root lines"

def als(free, r, iters=600, seed=0, ridge=1e-9):
    V = free.shape[0]; rng = np.random.default_rng(seed)
    X = rng.standard_normal((V, r)); Y = rng.standard_normal((V, r))
    req = ~free
    I = np.eye(r)
    def solve(Y):
        # row a: min sum_{b req} (x.Y_b - delta_ab)^2  -> (sum_{b req} Y_b Y_b^T) x = Y_a
        G = Y.T @ Y
        A = np.empty((V, r, r))
        for a in range(V):
            F = Y[free[a]]
            A[a] = G - F.T @ F + ridge * I
        return np.linalg.solve(A, Y[:, :, None])[:, :, 0]
    best = None
    for it in range(iters):
        X = solve(Y); Y = solve(X)
        if it % 25 == 24 or it == iters - 1:
            N = X @ Y.T
            res = np.sum(((N - np.eye(V)) * req) ** 2)
            best = res
            if res < 1e-18: break
            # rescale to balance
            s = np.sqrt(np.linalg.norm(Y, axis=0) / (np.linalg.norm(X, axis=0) + 1e-30))
            X *= s; Y /= s
    return best

def main():
    if sys.argv[1] == 'e8':
        free, theta, h, name = e8_family(); rest = sys.argv[2:]
    else:
        n, q = int(sys.argv[2]), int(sys.argv[3])
        free, theta, h, name = pg_family(n, q); rest = sys.argv[4:]
    V = free.shape[0]; deg = free.sum(1)
    print(f"{name}: V={V}, theta={theta:.4f}, free degree {deg.min()}..{deg.max()}, theta*V={theta*V:.1f}", flush=True)
    rmin, rmax = (int(rest[0]), int(rest[1])) if rest else (1, V)
    tries = int(rest[2]) if len(rest) > 2 else 4
    iters = int(rest[3]) if len(rest) > 3 else 600
    feas = {}
    lo, hi = rmin, rmax          # find smallest r with residual ~ 0 (bisection, assumes monotone)
    while lo < hi:
        r = (lo + hi) // 2
        res = min(als(free, r, iters=iters, seed=s) for s in range(tries))
        ok = res < 1e-12; feas[r] = res
        print(f"  r={r}: best residual {res:.3e} {'FEASIBLE' if ok else ''}", flush=True)
        if ok: hi = r
        else: lo = r + 1
    mr = lo
    print(f"  numerical real minrank ~ {mr}  (mr/V = {mr/V:.4f})")
    print(f"  sandwich: {mr/V:.4f} <= gamma_min <= {mr/(theta*V):.4f}   [field-line ports]")
    print(f"  F2 minrank floor theta*V = {theta*V:.1f}; Q/F2 gap ratio mr/(theta V) = {mr/(theta*V):.3f}")

if __name__ == '__main__':
    main()

"""Translation-equivariant WH Tanner codes on biaffine(q), q = 2^m (char 0, exact Fourier formula).
Point code at every point: spectrum T_P (labels = slope m). Line (m,b): labels x -> x (affine relabels per
slope class are absorbed into a per-class spectrum T_m). Fourier over the translation group Z_2^{2m}:
  S0 = sum_{alpha,gamma} dim{ f in A_P : supp f subset M_{alpha,gamma} },
  M_{alpha,gamma} = { m : alpha + beta_m(gamma) in T_m },  beta_m(gamma)_i = <gamma, m*e_i>.
Trace bound: dim <= |M||T_P|/q, so S0/E <= 1/4 + O(1/q) for half-rate codes.
usage: equiv.py q iters seed  (hill-climb over T_P and the q sets T_m; prints best S0/E)"""
import sys, random, json
import numpy as np, galois

def setup(q):
    F = galois.GF(q); m = q.bit_length() - 1
    H = np.array([[(-1) ** bin(s & u).count('1') for u in range(q)] for s in range(q)], dtype=float)
    beta = np.zeros((q, q), dtype=int)  # beta[mm, gamma]
    for mm in range(q):
        for g in range(q):
            b = 0
            for i in range(m):
                if bin(g & int(F(mm) * F(1 << i))).count('1') & 1: b |= 1 << i
            beta[mm, g] = b
    return H, beta

def S0(q, H, beta, TP, Tm):
    TP = sorted(TP); inT = np.zeros((q, q), dtype=bool)
    for mm in range(q): inT[mm, list(Tm[mm])] = True
    tot = 0
    for g in range(q):
        for a in range(q):
            Mc = [mm for mm in range(q) if not inT[mm, a ^ beta[mm, g]]]
            if not Mc: tot += len(TP); continue
            tot += len(TP) - np.linalg.matrix_rank(H[np.ix_(TP, Mc)])
    return int(tot)

def climb(q, iters, seed):
    rng = random.Random(seed); H, beta = setup(q); k = q // 2
    TP = set(rng.sample(range(q), k)); Tm = [set(rng.sample(range(q), k)) for _ in range(q)]
    cur = S0(q, H, beta, TP, Tm); best = cur
    for it in range(iters):
        which = rng.randrange(q + 1)
        S = TP if which == q else Tm[which]
        a = rng.choice(sorted(S)); b = rng.choice([x for x in range(q) if x not in S])
        S.remove(a); S.add(b)
        new = S0(q, H, beta, TP, Tm)
        if new >= cur or rng.random() < 0.02:
            cur = new
        else:
            S.remove(b); S.add(a)
        if cur > best:
            best = cur; bestcfg = (sorted(TP), [sorted(t) for t in Tm])
    return best, (best / q ** 3)

if __name__ == '__main__':
    q, iters, seed = map(int, sys.argv[1:4])
    H, beta = setup(q)
    k = q // 2
    sub = [s for s in range(q) if not s & 1]
    print(json.dumps(dict(q=q, check_sub_sub=S0(q, H, beta, sub, [sub] * q), E=q ** 3)), flush=True)
    b, r = climb(q, iters, seed)
    print(json.dumps(dict(q=q, seed=seed, iters=iters, bestS0=b, SE=round(r, 4), thr=round(1 / 6, 4))), flush=True)
    print('PASS_EXPRF_EQUIV' if r < 1 / 6 else 'FAIL_EXPRF_EQUIV')   # the equivariant biaffine optimum found stays below the 1/6 threshold

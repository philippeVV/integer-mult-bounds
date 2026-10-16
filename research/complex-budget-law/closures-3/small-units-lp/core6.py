"""h = 6 odd-model core: isometry group of the standard form on F2^n, types (q, w) odd/odd with w.q = 1,
closed columns (W, Q) and their orbits under G x <swap>.  Cached to core{n}.pkl.
usage (build): core6.py n"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys, time, pickle, os, itertools
import numpy as np

def dot(a, b): return bin(a & b).count('1') & 1
def basis_of(vs):
    B = {}
    for x in vs:
        while x:
            h = x.bit_length() - 1
            if h in B: x ^= B[h]
            else: B[h] = x; break
    return list(B.values())
def span_mask(vs):
    S = [0]
    for b in basis_of(vs): S = S + [s ^ b for s in S]
    return sum(1 << s for s in S)
def bits(mask): return [i for i in range(mask.bit_length()) if (mask >> i) & 1]

def build(n):
    N = 1 << n
    ODD = sum(1 << v for v in range(N) if bin(v).count('1') & 1)
    odd = bits(ODD)
    # --- group generators as permutations of F2^n (linear isometries of the standard dot product)
    def lin(images):
        return tuple(np.bitwise_xor.reduce([images[i] for i in range(n) if (v >> i) & 1], initial=0) if v else 0 for v in range(N))
    def lin2(images):
        out = []
        for v in range(N):
            r = 0
            for i in range(n):
                if (v >> i) & 1: r ^= images[i]
            out.append(r)
        return tuple(out)
    e = [1 << i for i in range(n)]
    gens = [lin2([e[1], e[0]] + e[2:]), lin2(e[1:] + e[:1])]
    for a in [0b11, (1 << 4) - 1, N - 1] if n >= 4 else [0b11]:
        if bin(a).count('1') % 2 == 0 and a < N:
            gens.append(lin2([x ^ (a if dot(x, a) else 0) for x in e]))
    for g in gens:
        assert all(dot(g[x], g[y]) == dot(x, y) for x in range(N) for y in range(N))
    # group order by BFS (perms of N points)
    G = {tuple(range(N))}; fr = list(G)
    while fr:
        nf = []
        for p in fr:
            for g in gens:
                q = tuple(g[p[i]] for i in range(N))
                if q not in G: G.add(q); nf.append(q)
        fr = nf
    order = len(G)
    # --- types
    types = [(q, w) for q in odd for w in odd if dot(w, q)]
    tidx = {t: i for i, t in enumerate(types)}
    # --- odd-spanned subspaces, perp, closure
    subs = set(); fr = [1]; subs.add(1)
    while fr:
        nf = []
        for S in fr:
            for v in odd:
                if not (S >> v) & 1:
                    Sn = S
                    for s in bits(S): Sn |= 1 << (s ^ v)
                    if Sn not in subs: subs.add(Sn); nf.append(Sn)
        fr = nf
    subs = sorted(subs)
    PERP = {}
    def perp(S):
        if S not in PERP:
            bs = basis_of(bits(S))
            PERP[S] = sum(1 << v for v in range(N) if all(not dot(v, b) for b in bs))
        return PERP[S]
    PL = {}
    def plus(X, Y):
        k = (X, Y)
        if k not in PL: PL[k] = span_mask(basis_of(bits(X)) + basis_of(bits(Y)))
        return PL[k]
    SP = {}
    def spanodd(mask):
        if mask not in SP: SP[mask] = span_mask(bits(mask))
        return SP[mask]
    closed = set()
    for U in subs:
        for V in subs:
            Wm, Qm = U & ODD, V & ODD
            while True:
                W2 = plus(spanodd(Wm), perp(spanodd(Qm))) & ODD
                Q2 = plus(spanodd(Qm), perp(spanodd(W2))) & ODD
                if (W2, Q2) == (Wm, Qm): break
                Wm, Qm = W2, Q2
            closed.add((Wm, Qm))
    cols = sorted(closed)
    # beta = dim span W - dim(span W cap span(Q)^perp)
    def dim(S): return bin(S).count('1').bit_length() - 1
    beta = []
    for Wm, Qm in cols:
        A = spanodd(Wm); C = spanodd(Qm)
        beta.append(dim(A) - dim(A & perp(C)))
    # --- orbits by union-find under gens + swap
    def pmask(g, mask):
        r = 0
        for b in bits(mask): r |= 1 << g[b]
        return r
    cidx = {c: i for i, c in enumerate(cols)}
    par = list(range(len(cols)))
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for i, (Wm, Qm) in enumerate(cols):
        imgs = [(pmask(g, Wm), pmask(g, Qm)) for g in gens] + [(Qm, Wm)]
        for c in imgs:
            j = cidx[c]; a, b = find(i), find(j)
            if a != b: par[a] = b
    corb = [find(i) for i in range(len(cols))]
    tpar = list(range(len(types)))
    def tfind(x):
        while tpar[x] != x: tpar[x] = tpar[tpar[x]]; x = tpar[x]
        return x
    for i, (q, w) in enumerate(types):
        for t in [(g[q], g[w]) for g in gens] + [(w, q)]:
            a, b = tfind(i), tfind(tidx[t])
            if a != b: tpar[a] = b
    torb = [tfind(i) for i in range(len(types))]
    return dict(n=n, N=N, ODD=ODD, odd=odd, gens=gens, order=order, types=types, cols=cols, beta=beta, corb=corb, torb=torb)

def load(n):
    return build(n)   # closures-3: rebuilt on every call (the scouting run cached it to core{n}.pkl, about 2 MB)

if __name__ == '__main__':
    n = int(sys.argv[1]); t0 = time.time()
    D = load(n)
    from collections import Counter
    print('n=%d |G|=%d types=%d closed cols=%d (%.0f s)' % (n, D['order'], len(D['types']), len(D['cols']), time.time() - t0))
    to = Counter(D['torb']); print('type orbits:', len(to), sorted(to.values()))
    co = Counter(D['corb']); print('col orbits:', len(co))
    b = Counter(D['beta'][i] for i in set(D['corb'])); print('col orbits by beta:', sorted(b.items()))

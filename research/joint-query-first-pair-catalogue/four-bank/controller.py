#!/usr/bin/env python3
"""Literal replay of the whole retained nine-stage QQ16 replica (maps (1),(3) of CONSERVED-QQ16.md,
ordinary bridges, terminal) in exact rationals on a finite noncommuting permutation-route model,
with the full (generous) observation set:
  * every selected register of every bank after EVERY event once that bank's opening calls are
    complete (bank 0 after F0,B1; bank k>=1 after its first call B_{2k+1}); the transient F0 output
    and a bank's untouched inputs are not observations, exactly as in the reviewed contract,
  * the five forward histories f0..f8,
  * completed-X descendants = final Y_U of every bank (named final Y after (X,Y)->(Y,-X)),
  * optionally the PR2174 pre-F0 records q0=B0 X_0U, q2=B2 X_0U.
Checks the per-pair erased families:
  pair 0 : h in ker B1 (reviewed witness),
  pair k>=1 (first call B_{b,c}, b=2k+1, c=b+1, preceded by F_a on pair k-1, a=2k):
      h in ker B_a, X_U=h, Y_U=0, X_T=-(1/11)P_b B_b h, X_A=P_{c,A}(B_c h - ell_{c,T} X_T), Y_S=-X_S.
usage: controller.py ARCHIVE [m] [seed] [mutation]"""
import os, random, sys
from fractions import Fraction as Q
sys.path.insert(0, __file__.rsplit('/', 1)[0] + '/..')
import pair_catalogue as pc

A, T, U = pc.A, pc.T, pc.U
ALL = list(range(112))
m = int(sys.argv[2]) if len(sys.argv) > 2 else 5
rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 20261014)
MUT = sys.argv[4] if len(sys.argv) > 4 else ''
cols = pc.query_columns(pc.load(sys.argv[1]))
D = {0: [list(range(m))] * 112}
for u in range(1, 9):
    D[u] = [rng.sample(range(m), m) for _ in range(112)]
Z = lambda: [Q(0)] * m
def apply(P, v):
    w = Z()
    for i, x in enumerate(v):
        w[P[i]] += x
    return w
def applyinv(P, v):
    return [v[P[i]] for i in range(m)]
def add(*vs):
    return [sum(z) for z in zip(*vs)]
def sc(c, v):
    return [c * x for x in v]
def query(u, f, ports):
    out = [Z() for _ in range(8)]
    for p in ports:
        g = apply(D[u][p], f[p])
        for k in range(8):
            if cols[p][k]:
                out[k] = add(out[k], sc(cols[p][k], g))
    return out
def inject(u, S, q):
    f = {}
    for p in S:
        v = Z()
        for k in range(8):
            if cols[p][k]:
                v = add(v, sc(Q(cols[p][k], 4), q[k]))
        f[p] = applyinv(D[u][p], v)
    return f
def kU(u, q):
    return {p: applyinv(D[u][p], [sum(Q(cols[p][k], 4) * q[k][i] for k in range(8)) for i in range(m)]) for p in U}
def adjoint(u, q):
    return {p: applyinv(D[u][p], [sum(cols[p][k] * q[k][i] for k in range(8)) for i in range(m)]) for p in U}
def forward(a, x, y):
    f = query(a, x, ALL)
    xa, xt, r = inject(a, A, query(a, x, T)), inject(a, T, f), kU(a, f)
    X, Y = {}, {}
    for p in U:
        X[p] = add(x[p], sc(-1, r[p])); Y[p] = add(y[p], X[p])
    for p in A + T:
        X[p] = xa[p] if p in A else xt[p]; Y[p] = add(x[p], y[p], sc(-1, X[p]))
    return X, Y, f
def backward(s, t, x, y):
    ua = inject(t, A, query(t, y, ALL))
    u = {p: (ua[p] if p in A else y[p]) for p in ALL}
    a, b = query(s, x, T), query(s, u, ALL)
    r = kU(s, a)
    xt = inject(s, T, [add(b[k], sc(12, a[k])) for k in range(8)])
    sub = inject(s, T, [add(b[k], sc(11, a[k])) for k in range(8)])
    X, Y = {}, {}
    for p in U:
        X[p] = add(x[p], sc(-1, y[p])); Y[p] = add(y[p], r[p])
    for p in A:
        X[p] = add(x[p], y[p], sc(-1, u[p])); Y[p] = u[p]
    for p in T:
        X[p] = xt[p]; Y[p] = add(u[p], sc(-1, sub[p]))
    return X, Y

# Event list: retained receipt (QQ16-SELECTED-WITNESS.json) with recipe pairs 0,0,0,1,1,2,2,3,3.
PAIR = [0, 0, 0, 1, 1, 2, 2, 3, 3] if MUT != 'pairmod' else [0, 0, 1, 1, 2, 2, 3, 3, 0]
EVENTS = [('local', 0), ('local', 1),
          ('bridge', 1, 0, 1), ('local', 2), ('bridge', 1, 0, -1), ('local', 3),
          ('bridge', 2, 1, 1), ('local', 4), ('bridge', 2, 1, -1), ('local', 5),
          ('bridge', 3, 2, 1), ('local', 6), ('bridge', 3, 2, -1), ('local', 7), ('local', 8),
          ('bridge', 2, 3, 1), ('bridge', 1, 2, 1), ('bridge', 0, 1, 1), ('terminal',)]
def recipe_events():
    import tarfile, gzip, json, hashlib
    raw = open(sys.argv[1], 'rb').read()
    pc.need(hashlib.sha256(raw).hexdigest() == pc.ARCHIVE_SHA, 'archive pin')
    import io
    with tarfile.open(fileobj=io.BytesIO(raw), mode='r:gz') as ar:
        rec = json.loads(gzip.decompress(ar.extractfile('t16/global/base/split-centres/GLOBAL-TEMPLATE.json.gz').read()))['recipe']
    ev, pairs = [], {}
    for e in rec:
        if e['kind'] == 'local':
            pc.need(e['reverse'] == (e['stage'] % 2 == 1), 'B on odd stages'); ev.append(('local', e['stage'])); pairs[e['stage']] = e['pair']
        elif e['kind'] == 'bridge':
            ev.append(('bridge', e['src'], e['dst'], e['coefficient']))
        else:
            pc.need(e == {'kind': 'terminal', 'pairs': 4}, 'terminal'); ev.append(('terminal',))
    return ev, [pairs[s] for s in range(9)]
_ev, _pairs = recipe_events()
pc.need(_ev == EVENTS and _pairs == [0, 0, 0, 1, 1, 2, 2, 3, 3], 'event list bound to the pinned global recipe')
if MUT == 'nobridge':
    EVENTS = [e for e in EVENTS if e != ('bridge', 1, 0, -1)]

def run(state):
    """state: list of 4 banks (X,Y). Returns observation vector (list of Fractions) and final state."""
    X = [dict(b[0]) for b in state]; Y = [dict(b[1]) for b in state]
    obs = []
    done = set()   # banks whose opening call(s) are complete: bank 0 after B1, bank k after B_{2k+1}
    FINAL = os.environ.get('QQ16_CONTRACT') in ('final', 'finalf')
    KEEPF = os.environ.get('QQ16_CONTRACT') != 'final'
    def cut():
        if FINAL:
            return
        for k in sorted(done):
            for p in A + T:
                obs.extend(X[k][p]); obs.extend(Y[k][p])
    for ev in EVENTS:
        if ev[0] == 'local':
            s = ev[1]; k = PAIR[s]
            if s % 2 == 0:
                X[k], Y[k], f = forward(s, X[k], Y[k])
                if KEEPF:
                    for r in f:
                        obs.extend(r)
            else:
                X[k], Y[k] = backward(s, s + 1, X[k], Y[k])
                done.add(k)
            cut()
        elif ev[0] == 'bridge':
            _, src, dst, c = ev
            for p in U:
                X[dst][p] = add(X[dst][p], sc(c, X[src][p]))
            for p in U:
                Y[src][p] = add(Y[src][p], sc(-c, Y[dst][p]))
            cut()
        else:
            for k in range(4):
                X[k], Y[k] = Y[k], {p: sc(-1, v) for p, v in X[k].items()}
            cut()
    if FINAL:
        for k in range(4):
            for p in A + T:
                obs.extend(X[k][p]); obs.extend(Y[k][p])
    for k in range(4):
        for p in U:
            obs.extend(Y[k][p])          # completed-X descendants (named final Y)
    return obs, X, Y

def zero_bank():
    return ({p: Z() for p in ALL}, {p: Z() for p in ALL})

def nz(v):
    return any(x != 0 for x in v)

def kernel_vec(a):
    v = {p: [Q(rng.randint(-3, 3)) for _ in range(m)] for p in U}
    bv = adjoint(a, query(a, v, U))
    return {p: add(sc(48, v[p]), sc(-1, bv[p])) for p in U}

def witness(k):
    banks = [zero_bank() for _ in range(4)]
    x, y = banks[k]
    if k == 0:
        h = kernel_vec(1)
        assert not any(nz(r) for r in query(1, h, U))
        B0h, B2h = query(0, h, U), query(2, h, U)
        i0, i2 = inject(0, A, B0h), inject(2, A, B2h)
        xt = inject(0, T, query(0, i2, A))
        for p in U: x[p] = h[p]
        for p in A: x[p] = sc(-1, add(i0[p], i2[p])); y[p] = sc(-1, x[p])
        for p in T: x[p] = xt[p]; y[p] = sc(-1, xt[p])
    else:
        a, b, c = 2 * k, 2 * k + 1, 2 * k + 2
        h = kernel_vec(a)
        if MUT == 'notker':
            h = {p: [Q(rng.randint(-3, 3)) for _ in range(m)] for p in U}
        assert MUT == 'notker' or not any(nz(r) for r in query(a, h, U))
        fac = Q(-1, 11) if MUT != 'no11' else Q(-1, 12)
        xt = inject(b, T, [sc(fac, r) for r in query(b, h, U)])
        lt = query(c, xt, T)
        xa = inject(c, A, [add(r1, sc(-1, r2)) for r1, r2 in zip(query(c, h, U), lt)])
        for p in U: x[p] = h[p]
        for p in A: x[p] = xa[p]; y[p] = sc(-1, xa[p])
        for p in T: x[p] = xt[p]; y[p] = sc(-1, xt[p])
    J = [v for p in A + T for v in x[p]]
    q0 = query(0, x if k == 0 else banks[0][0], U); q2 = query(2, banks[0][0], U)
    return banks, J, q0, q2

def main():
  total = None
  for k in range(4):
    banks, J, q0, q2 = witness(k)
    obs, X, Y = run(banks)
    zero = not nz(obs)
    print(f'pair {k}: selected input nonzero={nz(J)}  all observations zero={zero}  '
          f'PR2174 pair-0 records zero={not nz([v for r in q0 + q2 for v in r])}  hole nonzero={any(nz(X[j][p]) for j in range(4) for p in U)}')
    assert nz(J) and zero, 'witness failed'
    total = banks if total is None else [tuple({p: add(tb[i][p], b[i][p]) for p in ALL} for i in range(2)) for tb, b in zip(total, banks)]
  obs, _, _ = run(total)
  assert not nz(obs)
  print('PASS_QQ16_FOUR_PAIR_ERASURE m=%d sum-of-four zero=True' % m)

import os
if not os.environ.get('QQ16_LIB'):
    main()

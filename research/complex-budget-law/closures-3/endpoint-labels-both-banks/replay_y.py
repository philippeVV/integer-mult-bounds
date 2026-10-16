#!/usr/bin/env python3
"""expQW replay_y.py: expNP's gate-level replay extended to Y labels (both banks), one block h = 9, E8 P112.

Labels: the 9-bit E8 line labels (standard dot form; p_a.p_b = 1 + (r_a.r_b) mod 2), P112 = P120 minus an
orthogonal frame, sigma = an orthogonal perfect matching, q_b = p_sigma(b).
Endpoints: X_b source <p_b>, sink q_b^perp; Y_b source <q_b>, sink p_b^perp; scratch source 0, sink F.
Exact Fraction values over inputs x_b, y_b, z_r. Every gate: common nondegenerate frame, nested steps
(the engine routes a non-nested move as drop-to-a-nondegenerate-subspace then climb, and charges the drop).
Checks per schedule: nondegenerate frames, nested steps, exact (x, y, z) -> (-y, x, z), deficit identity
s = Wm - 2 sum u + 2L, data-drop bookkeeping, loss per window.
Schedules:
  Y0  data-drop baseline: L = 2V = sum u exactly (delta = 0)
  Y1  private helpers in all three windows (no data drops)
  (A dense-w1 schedule was drafted but dropped: its scratch z-cancellation was not finished in the time box.
   Lemma I does not need it -- flat.py certifies the span fact it rests on.)
Run: $S/expO/run replay_y.py  (needs networkx; ~1-2 min)
"""
import sys, os
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import networkx as nx
from e8 import frames, LAB, dotq, NL, IP

OK = True
def check(c, msg):
    global OK
    print(('PASS ' if c else 'FAIL ') + msg); OK &= bool(c)

H = 9
def basis(vecs):
    piv = {}
    for r in vecs:
        while r:
            h = r.bit_length() - 1
            if h in piv: r ^= piv[h]
            else: piv[h] = r; break
    keys = sorted(piv)
    for k in keys:
        for k2 in keys:
            if k2 != k and piv[k2] >> k & 1: piv[k2] ^= piv[k]
    return tuple(sorted(piv.values()))
def span_all(U):
    S = [0]
    for u in U: S += [s ^ u for s in S]
    return S
def contains(U, v): return basis(list(U) + [v]) == U
def subset(U, W): return all(contains(W, u) for u in U)
def perp(U): return basis([v for v in range(1, 1 << H) if all(dotq(v, u) == 0 for u in U)])
def nondeg(U): return not any(all(dotq(u, w) == 0 for w in U) for u in span_all(U) if u)
def meet(U, W): return basis([v for v in span_all(U) if contains(W, v)])
FULL = basis([1 << i for i in range(H)]); ZERO = ()

F0 = set(frames()[0]); P = [a for a in range(NL) if a not in F0]; V = len(P)
G = nx.Graph(); G.add_nodes_from(P)
G.add_edges_from((a, b) for i, a in enumerate(P) for b in P[i + 1:] if dotq(LAB[a], LAB[b]) == 0)
sig = {}
for a, b in nx.max_weight_matching(G, maxcardinality=True): sig[a] = b; sig[b] = a
check(len(sig) == V and all(dotq(LAB[a], LAB[sig[a]]) == 0 for a in P), f'P112: orthogonal perfect matching sigma ({V // 2} pairs)')
p = {a: LAB[a] for a in P}; q = {a: LAB[sig[a]] for a in P}
PL = {a: basis([p[a]]) for a in P}; QL = {a: basis([q[a]]) for a in P}
PP = {a: perp(PL[a]) for a in P}; QP = {a: perp(QL[a]) for a in P}
check(all(nondeg(x) for d in (PL, QL, PP, QP) for x in d.values()), 'all endpoint frames nondegenerate')

def best_sub(f, g):
    """largest nondegenerate subspace of f cap g that contains... (greedy: meet if nondeg, else search down)."""
    m = meet(f, g)
    if nondeg(m): return m
    best = ZERO
    for v in span_all(m):           # nondegenerate hyperplanes / lines of m: try v^perp cap m and <v>
        if not v: continue
        for cand in (meet(m, perp(basis([v]))), basis([v])):
            if nondeg(cand) and len(cand) > len(best): best = cand
    return best

class Net:
    def __init__(s, name):
        s.name, s.frame, s.val, s.src, s.kind = name, {}, {}, {}, {}
        s.s = s.L = 0; s.win = {}; s.cur = None; s.bad = []; s.ddrops = 0
        for b in P:
            s.new(('X', b), PL[b], {('x', b): Fr(1)}, 'X'); s.new(('Y', b), QL[b], {('y', b): Fr(1)}, 'Y')
    def new(s, r, fr, val, kind): s.frame[r], s.val[r], s.src[r], s.kind[r] = fr, dict(val), fr, kind
    def scratch(s, r): s.new(r, ZERO, {('z',) + r: Fr(1)}, 'S')
    def step(s, r, g):
        f = s.frame[r]
        if not nondeg(g): s.bad.append(('degenerate', r))
        if subset(f, g): s.s += len(g) - len(f)
        elif subset(g, f):
            d = len(f) - len(g); s.s += d; s.L += d; s.win[s.cur] = s.win.get(s.cur, 0) + d
            if s.kind[r] != 'S': s.ddrops += d; s.bad.append(('data drop', r))
        else: s.bad.append(('non-nested', r))
        s.frame[r] = g
    def goto(s, r, g):                       # route: drop to a nondegenerate subspace of f cap g, then climb
        f = s.frame[r]
        if not subset(f, g) and not subset(g, f): s.step(r, best_sub(f, g))
        s.step(r, g)
    def gate(s, g, regs, ops):
        for r in regs: s.goto(r, g)
        for t, sr, c in ops:
            for k, v in list(s.val[sr].items()):
                s.val[t][k] = s.val[t].get(k, 0) + c * v
                if s.val[t][k] == 0: del s.val[t][k]
    def join(s, regs, extra=ZERO):           # common frame = sum of current frames (+extra); nondeg check by step
        return basis([v for r in regs for v in s.frame[r]] + list(extra))
    def finish(s, allow_data_drops=False):
        s.cur = 'sink'
        for r in list(s.frame):
            k = s.kind[r]; sink = FULL if k == 'S' else (QP[r[1]] if k == 'X' else PP[r[1]])
            s.goto(r, sink)
        Wn = len(s.frame); su = sum(len(s.src[r]) for r in s.frame)
        wrong = [r for r in s.frame if s.val[r] != ({('y', r[1]): Fr(-1)} if s.kind[r] == 'X' else {('x', r[1]): Fr(1)} if s.kind[r] == 'Y' else {('z',) + r: Fr(1)})]
        if wrong: print('   wrong values:', len(wrong), [(r, s.val[r]) for r in wrong[:2]])
        okv = all(s.val[('X', b)] == {('y', b): Fr(-1)} and s.val[('Y', b)] == {('x', b): Fr(1)} for b in P) and \
            all(s.val[r] == {('z',) + r: Fr(1)} for r in s.frame if s.kind[r] == 'S')
        bad = [x for x in s.bad if not (allow_data_drops and x[0] == 'data drop')]
        check(not bad, f'{s.name}: frames nondegenerate, steps nested' + ('' if allow_data_drops else ', no data drop') + (f' {bad[:3]}' if bad else ''))
        check(okv, f'{s.name}: exact (x, y, z) -> (-y, x, z), scratch restored')
        check(s.s == Wn * H - 2 * su + 2 * s.L, f'{s.name}: deficit identity {s.s} = {Wn}*9 - 2*{su} + 2*{s.L}')
        D = Wn * H - s.s
        print(f'   {s.name}: W = {Wn}, sum u = {su}, L = {s.L} (data drops {s.ddrops}), windows {s.win}, D = {D}, delta = {D / (Wn * H):+.5f}')
        return s.L, su, D

def private(net, tag, xy):
    """one private helper per port; xy: y_b += x_b, else x_b -= y_b. Helper bounces between <p_b> and <q_b> frames."""
    net.cur = tag
    for b in P:
        c = (tag, 'c', b); net.scratch(c); X, Y = ('X', b), ('Y', b)
        if xy:
            net.gate(net.frame[Y], [Y, c], [(Y, c, -1)]); net.gate(net.frame[X], [X, c], [(c, X, 1)])
            net.gate(net.frame[Y], [Y, c], [(Y, c, 1)]); net.gate(net.frame[X], [X, c], [(c, X, -1)])
        else:
            net.gate(net.frame[X], [X, c], [(X, c, 1)]); net.gate(net.frame[Y], [Y, c], [(c, Y, 1)])
            net.gate(net.frame[X], [X, c], [(X, c, -1)]); net.gate(net.frame[Y], [Y, c], [(c, Y, -1)])

def datadrop(net):
    net.cur = 'dd'
    for b in P:
        X, Y = ('X', b), ('Y', b)
        net.gate(PL[b], [X, Y], [(Y, X, 1), (X, Y, -1), (Y, X, 1)])   # Y drops <q_b> -> 0 and climbs to <p_b>

# ---- Gram decoder on P112: N = G/2 = sum_k c_k c_k^T / 4, c_b(k) = r_b . f_k  (f_k the removed frame)
FR = sorted(F0)
c = {b: [IP[b][f] for f in FR] for b in P}
N = {(a, b): Fr(sum(c[a][k] * c[b][k] for k in range(8)), 4) for a in P for b in P}
check(all(N[a, a] == 1 for a in P) and all(N[a, b] == 0 for a in P for b in P if a != b and dotq(p[a], p[b]) == 1),
      'N = G/2 on P112 is a masked identity (8 terms)')
S_k = [[b for b in P if c[b][k]] for k in range(8)]
CORR = [(a, b) for a in P for b in P if a != b and N[a, b] != 0]
check(all(dotq(p[a], p[b]) == 0 for a, b in CORR), f'{len(CORR)} corrections, all on free (orthogonal) pairs')
print(f'   term spans: {[len(basis([p[b] for b in S])) for S in S_k]}, pairing ranks: '
      f'{[len(basis([sum(dotq(p[a], p[b]) << i for i, b in enumerate(S)) for a in S])) for S in S_k]}')

print('\n== Y0 data-drop baseline'); n0 = Net('Y0'); datadrop(n0); R0 = n0.finish(allow_data_drops=True)
check(R0[0] == R0[1] and R0[2] == 0, 'Y0: L = sum u = 2V, D = 0')
print('\n== Y1 private helpers, three windows'); n1 = Net('Y1')
private(n1, 'w1', True); private(n1, 'w2', False); private(n1, 'w3', True); R1 = n1.finish()
print('\n== summary (V = %d, sum u = %d; a witness needs L < sum u)' % (V, 2 * V))
for nm, R in (('Y0', R0), ('Y1', R1)): print(f'   {nm}: L = {R[0]}, D = {R[2]}')
check(min(R0[0], R1[0]) >= 2 * V, 'no replayed both-bank schedule reaches L < sum u')
print('\nALL PASS' if OK else '\nSOME CHECK FAILED'); print('PASS_EXPQW_REPLAY' if OK else 'FAIL_EXPQW_REPLAY'); sys.exit(0 if OK else 1)

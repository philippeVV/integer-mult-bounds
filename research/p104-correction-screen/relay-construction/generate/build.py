"""Build an explicit X-only correction word from a pool design (requests child->parent with forms).
Ops (all logical, exact Fractions):
  ['move', L, frame_basis]          lifetime L transports to a strictly larger frame (or ZERO->first frame)
  ['read', L, s, c, cut]            L += c * X_s  (X_s frame must equal L's frame at that cut)
  ['shear', i, j, c]                L_i += c * L_j (equal frames)
  ['write', t, L, c, cut]           Y_t += c * L  (Y_t frame must equal L's frame at that cut)
Phases: cut 0 (compensation at ZERO, reads on lines, pool mixing), cut NE (arrivals, writes, FULL cleanup)."""
import json, sys, collections
from fractions import Fraction as Q
geo = json.load(open('geo.json')); n = geo['n']; q = geo['q']; NE = geo['NE']
des = json.load(open(sys.argv[1]))
def basis_of(m):
    b = []
    for v in range(512):
        if m >> v & 1:
            w = v
            for x in b: w = min(w, w ^ x)
            if w: b.append(w)
    return b
dimf = lambda m: (m.bit_count()).bit_length() - 1
req = [(int(ch), (tuple(par) if isinstance(par, list) else int(par)), {int(k): Q(v) for k, v in f.items()}) for ch, par, f in des['req']]
pools = sorted({ch for ch, _, _ in req} | {p for _, p, _ in req if not isinstance(p, tuple)}, key=lambda m: (dimf(m), m))
inc = collections.defaultdict(list); out = collections.defaultdict(list)
for k, (ch, par, f) in enumerate(req):
    out[ch].append(k)
    if not isinstance(par, tuple): inc[par].append(k)
ops = []; content = {}; frame = {}; nlife = 0
def new_life():
    global nlife
    L = nlife; nlife += 1; content[L] = {('d', L): Q(1)}; frame[L] = 1; return L
def addto(a, c, b):  # content[a] += c * content[b]
    for k, v in content[b].items():
        z = content[a].get(k, 0) + c*v
        if z: content[a][k] = z
        else: content[a].pop(k, None)
def xpart(L): return {k[1]: v for k, v in content[L].items() if k[0] == 'x'}
def solve(G, f):
    """coefficients a with sum a_i G_i == f (G list of dicts); exact Gaussian elimination"""
    keys = sorted({k for g in G for k in g} | set(f)); m = len(G)
    A = [[g.get(k, Q(0)) for g in G] + [f.get(k, Q(0))] for k in keys]
    piv = []; r = 0
    for c in range(m):
        p = next((i for i in range(r, len(A)) if A[i][c]), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]; A[r] = [x / A[r][c] for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c]: A[i] = [x - A[i][c]*y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    if any(row[-1] for row in A[r:]): raise SystemExit('BUILD FAIL: form not in span')
    a = [Q(0)]*m
    for i, c in enumerate(piv): a[c] = A[i][-1]
    return a
scale = {}
carrier = {}   # request index -> lifetime carrying its form (at child frame)
plan = []      # deferred op list in execution order (cut-0 part)
def mv(L, m):
    assert dimf(m) > dimf(frame[L]) and (frame[L] & m) == frame[L]
    frame[L] = m; plan.append(['move', L, basis_of(m)])
for V in pools:
    ents = [carrier[k] for k in inc[V]]
    for L in ents: mv(L, V)
    exits = out[V]
    if dimf(V) == 1:
        s = q.index(basis_of(V)[0])
        for k in exits:
            L = new_life(); mv(L, V); c = req[k][2][s]
            content[L][('x', s)] = c; plan.append(['read', L, s, c, 0]); carrier[k] = L; scale[k] = Q(1)
        continue
    members = list(ents); assigned = set()
    def indep(order):
        bs = []
        for L in order:
            g = xpart(L)
            if not g: continue
            if bs:
                try: solve([xpart(B) for B in bs], g); continue
                except SystemExit: pass
            bs.append(L)
        return bs
    for k in exits:
        f = {s_: v for s_, v in req[k][2].items() if v}
        order = [L for L in members if L not in assigned] + [L for L in members if L in assigned]
        bs = indep(order)
        a = solve([xpart(B) for B in bs], f)
        pick = next((B for B, c in zip(bs, a) if c and B not in assigned), None)
        if pick is None:
            L = new_life(); mv(L, V); members.append(L)
            for B, c in zip(bs, a):
                if c: addto(L, c, B); plan.append(['shear', L, B, c])
            scale[k] = Q(1)
        else:
            L = pick; cp = dict(zip(bs, a))[pick]
            for B, c in zip(bs, a):
                if c and B != pick: addto(L, c / cp, B); plan.append(['shear', L, B, c / cp])
            scale[k] = cp
        assert {s_: v * scale[k] for s_, v in xpart(L).items()} == f
        assigned.add(L); carrier[k] = L
# arrivals and writes at cut NE
Hframe = {t: geo['Y'][t][-1] for t in range(n)}
final_writes = []
FULL = (1 << 512) - 1
for k, (ch, par, f) in enumerate(req):
    if isinstance(par, tuple):
        t = par[1]; L = carrier[k]
        Hm = 0
        for v in range(512):
            pass
        final_writes.append((t, L, scale[k]))
# dirt compensation: Y_t receives content[L] for each arrival; subtract the d-part at cut 0 (all at ZERO)
comp = []
for t, L, sc in final_writes:
    for key, v in content[L].items():
        if key[0] == 'd': comp.append(['write', t, key[1], -v * sc, 0])
ops = [['move_zero_all', nlife]] + comp + plan
# cut NE: move arrivals to H_t, write, then all to FULL
def span_mask(bs):
    s = {0}
    for b in bs:
        if b not in s: s |= {x ^ b for x in s}
    m = 0
    for v in s: m |= 1 << v
    return m
for t, L, sc in final_writes:
    Hm = span_mask(Hframe[t]['basis'])
    if frame[L] != Hm: mv(L, Hm); ops.append(plan.pop())
    ops.append(['write', t, L, sc, NE])
for L in range(nlife):
    if frame[L] != FULL: mv(L, FULL); ops.append(plan.pop())
# restoration at FULL: undo shears in reverse, then subtract reads
shears = [o for o in ops if o[0] == 'shear']; reads = [o for o in ops if o[0] == 'read']
for o in reversed(shears): ops.append(['shear', o[1], o[2], -o[3]])
for o in reads: ops.append(['read', o[1], o[2], -o[3], NE])
prof = collections.Counter()
json.dump(dict(lifetimes=nlife, ops=[[str(x) if isinstance(x, Q) else x for x in o] for o in ops]), open(sys.argv[2], 'w'))
print('lifetimes', nlife, 'ops', len(ops), 'shears', len(shears), 'reads', len(reads), 'writes', sum(1 for o in ops if o[0] == 'write'))

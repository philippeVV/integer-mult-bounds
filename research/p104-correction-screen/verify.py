"""expAZ: structure, frozen-frame floors and budget for the P104 I-G/2 correction.
Standard library only.  Reads (read-only) files of this repository.
Usage: python3 -I -B verify_expAZ.py [--skip-orbits]   -> writes RESULTS.json next to this file."""
import json, sys, os, collections, itertools
from fractions import Fraction as Q
from itertools import combinations, product
HERE = os.path.dirname(os.path.abspath(__file__)); SRC = os.path.join(HERE, '..', '..')
sys.path.insert(0, os.path.join(SRC, 'research', 'complex-local-word-ceiling'))
import ceiling as CEIL
R = {}
def need(ok, why):
    if not ok: raise SystemExit('FAIL: ' + why)
# ---------------- A. algebraic structure ----------------
def e8_lines():
    L = set()
    for i, j in combinations(range(8), 2):
        for sg in (-1, 1):
            v = [0]*8; v[i] = 2; v[j] = 2*sg; L.add(tuple(v))
    for sg in product((-1, 1), repeat=7):
        if sg.count(-1) % 2 == 0: L.add((1, *sg))
    return sorted(L)
rr = e8_lines(); need(len(rr) == 120, '120 lines')
fx = json.load(open(os.path.join(SRC, 'research/p104-native-channel/inputs/P104-FIXTURE.json')))['section4_explicit_P104']
removed = {tuple(r) for r in fx['F0'] + fx['C']}
ports = [r for r in rr if r not in removed]; n = len(ports); need(n == 104, 'P104')
ip = lambda a, b: sum(x*y for x, y in zip(a, b)) // 4
G = [[ip(a, b) for b in ports] for a in ports]
mm = lambda A, B: [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
d = lambda i, j: int(i == j)
G2 = mm(G, G); need(all(G2[i][j] == 26*G[i][j] for i in range(n) for j in range(n)), 'G^2=26G')
S = [[G[i][j] - 2*d(i, j) for j in range(n)] for i in range(n)]
S2 = mm(S, S); need(all(S2[i][j] == 22*S[i][j] + 48*d(i, j) for i in range(n) for j in range(n)), 'S^2=22S+48I')
C = [[Q(-S[i][j], 2) for j in range(n)] for i in range(n)]          # C = I - G/2 = -S/2
need(all(C[i][j] == d(i, j) - Q(G[i][j], 2) for i in range(n) for j in range(n)), 'C=I-G/2')
Cinv = [[d(i, j) - Q(G[i][j], 24) for j in range(n)] for i in range(n)]
need(mm(C, Cinv) == [[Q(d(i, j)) for j in range(n)] for i in range(n)], 'inverse I-G/24')
def rank_mod(M, p=None):
    M = [[(Q(x) if p is None else x % p) for x in r] for r in M]; rk = 0
    for c in range(len(M[0])):
        pv = next((i for i in range(rk, len(M)) if M[i][c]), None)
        if pv is None: continue
        M[rk], M[pv] = M[pv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = (M[i][c]/M[rk][c]) if p is None else M[i][c]*pow(M[rk][c], -1, p) % p
                M[i] = [(a - f*b) if p is None else (a - f*b) % p for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk
Pm = [[Q(G[i][j], 26) for j in range(n)] for i in range(n)]
resid = [[C[i][j] + 12*Pm[i][j] for j in range(n)] for i in range(n)]      # C - E, E = -12P (rank 8)
need(resid == [[d(i, j) - Pm[i][j] for j in range(n)] for i in range(n)], 'C+12P = I-P')
R['A_structure'] = dict(
    G_sq='26G', S_sq='22S+48I', C='-S/2 = (I-P) - 12P, P=G/26 rank-8 projector',
    C_inverse='I - G/24', eigen_C={'1': 96, '-12': 8}, det_C='(-12)^8',
    rank_Q_C=rank_mod(C), rank_Q_residual_I_minus_P=rank_mod(resid), rank_Q_P=rank_mod(Pm),
    support_C=sum(1 for i in range(n) for j in range(n) if C[i][j]), row_support=48,
    ranks_S={'F2': rank_mod(S, 2), 'F3': rank_mod(S, 3), 'F5': rank_mod(S, 5), 'F13': rank_mod(S, 13)})
# Smith normal form of S via determinantal-divisor-free elimination
def snf(M):
    M = [r[:] for r in M]; rows, cols = len(M), len(M[0]); out = []; t = 0
    while t < min(rows, cols):
        nz = [(abs(M[i][j]), i, j) for i in range(t, rows) for j in range(t, cols) if M[i][j]]
        if not nz: break
        _, i, j = min(nz); M[t], M[i] = M[i], M[t]
        for r in M: r[t], r[j] = r[j], r[t]
        while True:
            for i in range(t+1, rows):
                qq = M[i][t] // M[t][t]
                if qq: M[i] = [a - qq*b for a, b in zip(M[i], M[t])]
            for j in range(t+1, cols):
                qq = M[t][j] // M[t][t]
                if qq:
                    for r in M: r[j] -= qq*r[t]
            rest = [(abs(M[i][j]), i, j) for i in range(t, rows) for j in range(t, cols) if (i == t or j == t) and (i, j) != (t, t) and M[i][j]]
            if rest:
                _, i, j = min(rest); M[t], M[i] = M[i], M[t]
                for r in M: r[t], r[j] = r[j], r[t]
                continue
            bad = [(i, j) for i in range(t+1, rows) for j in range(t+1, cols) if M[i][j] % M[t][t]]
            if not bad: break
            M[t] = [a + b for a, b in zip(M[t], M[bad[0][0]])]
        out.append(abs(M[t][t])); t += 1
    return dict(sorted(collections.Counter(out).items()))
R['A_structure']['smith_S'] = snf(S); R['A_structure']['smith_G'] = snf(G)
need(R['A_structure']['smith_S'] == {1: 8, 2: 88, 48: 8}, 'SNF')
# B. masked line-port model: supp C is exactly the B=0 off-diagonal set; diagonal has B=1
B = [[(G[i][j] + 1) % 2 for j in range(n)] for i in range(n)]
need(all((C[i][j] != 0) == (B[i][j] == 0 and i != j) for i in range(n) for j in range(n)), 'supp C = {B=0}')
need(all(B[i][i] == 1 for i in range(n)), 'B diagonal')
R['B_masked_model'] = 'supp(C) = {(t,s): t!=s, B_ts=0}: every term of C is a zero-pairing-rank (rank-zero cross-addition) term; the diagonal (B=1) is the decoder G/2, so the whole I-G/2 cost is chronological.'
# ---------------- C. frozen-frame chronology of the native channel ----------------
W = json.load(open(os.path.join(SRC, 'research/p104-native-channel/inputs/LOWRANK-WITNESS.json')))
def canon(vs):
    basis = []
    for v in vs:
        for b in basis: v = min(v, v ^ b)
        if v: basis = [min(b, b ^ v) for b in basis] + [v]
    return tuple(sorted(basis))
def span(fr):
    s = {0}
    for b in fr: s |= {x ^ b for x in s}
    return frozenset(s)
cur = {r: canon(f) for r, f in W['initial_frames'].items()}; hist = collections.defaultdict(list)
for r, f in cur.items(): hist[r].append((0, f))
for i, e in enumerate(W['program']):
    if e[0] == 'zero': cur['C'] = canon(e[2]); continue
    if e[0] in ('copy', 'erase'): continue
    if e[0] == 'transport':
        need(cur[e[1]] == canon(e[2]), 'transport entry'); cur[e[1]] = canon(e[3]); hist[e[1]].append((i+1, cur[e[1]]))
    else: need(cur[e[1]] == cur[e[2]] == canon(e[4]), 'common chart add')
X = [f'X{i}' for i in range(n)]; Y = [f'Y{i}' for i in range(n)]
for r in X + Y: need(cur[r] == canon(W['final_frames'][r]), 'final frame')
lastX = max(h[-1][0] for x in X for h in [hist[x]]); firstY = min(hist[y][1][0] for y in Y)
q = {s: hist[X[s]][0][1][0] for s in range(n)}
SF = collections.defaultdict(set); TF = collections.defaultdict(set)
for s in range(n):
    for _, f in hist[X[s]]: SF[f].add(s)
for t in range(n):
    for _, f in hist[Y[t]]:
        if f: TF[f].add(t)
SS = {v: span(v) for v in SF}; TS = {w: span(w) for w in TF}
Hs = {t: TS[hist[Y[t]][-1][1]] for t in range(n)}
inH = {(t, s) for t in range(n) for s in range(n) if q[s] in Hs[t]}
# port order of the witness equals the fixture port order (checked by the native package); bind Gamma:
need(inH == {(t, s) for t in range(n) for s in range(n) if C[t][s] != 0}, 'q_s in H_t  <=>  C_ts != 0')
nested = [(v, w) for v in SF for w in TF if SS[v] <= TS[w]]
cov = {(t, s) for v, w in nested for t in TF[w] for s in SF[v]}
need(cov == inH, 'nested access covers exactly supp C')
down = {v: [u for u in SF if SS[u] < SS[v]] for v in SF}; up = {w: [u for u in TF if TS[w] < TS[u]] for w in TF}
def unions(start, nbr, sets):
    out = set()
    def rec(f, acc):
        out.add(acc)
        for g in nbr[f]: rec(g, acc | sets[g])
    rec(start, frozenset(sets[start])); return out
Au = {v: unions(v, down, SF) for v in SF}; Bu = {w: unions(w, up, TF) for w in TF}
mcap = collections.defaultdict(int); best = 0
for v, w in nested:
    for a in Au[v]:
        for b in Bu[w]:
            P = [(t, s) for t in b for s in a if (t, s) in inH]
            best = max(best, len(P))
            for p in P: mcap[p] = max(mcap[p], len(P))
frac = sum(Q(1, mcap[p]) for p in inH)
# dimension-8 relay floor
lab = [q[s] for s in range(n)]; tot8 = 0; per8 = collections.Counter()
for t in range(n):
    ch = [f for _, f in hist[Y[t]] if f]; H = span(ch[-1]); Nt = {s for s in range(n) if lab[s] in H}
    low = set()
    for f in ch[:-1]:
        sp = span(f); low |= {s for s in Nt if lab[s] in sp}
    atH = {s for v, Ss in SF.items() if SS[v] == H for s in Ss}
    needs = Nt - low - atH; basis = list(ch[-1]); covers = []
    for mask in range(1, 256):
        Yset = set()
        for c in range(256):
            v = 0; par = 0
            for i in range(8):
                if c >> i & 1: v ^= basis[i]; par ^= (mask >> i) & 1
            if par == 0: Yset.add(v)
        covers.append({s for s in needs if lab[s] in Yset})
    k = 1
    while not any(set().union(*cc) >= needs for cc in combinations(covers, k)): k += 1
    tot8 += k; per8[k] += 1
R['C_frozen_frames'] = dict(
    last_source_transport_event=lastX, first_target_transport_event=firstY,
    source_frames=len(SF), target_frames=len(TF), nested_source_in_target_frame_pairs=len(nested),
    nested_access_covers_exactly_supp_C=True,
    no_relay_max_pairs_per_carrier=best,
    pairs_coverable_only_by_singleton_carriers=sum(1 for p in inH if mcap[p] == 1),
    no_relay_fractional_floor=str(frac), no_relay_carrier_floor=-(-frac.numerator // frac.denominator),
    relay_floor_distinct_dim8_visits=tot8, relay_floor_per_target=dict(per8))
# ---------------- D. budget (Colkitt ceiling engine, exact intervals) ----------------
kap = Q(836330098538459, 500000000000000000); asel = kap/(1-kap); abit = Q(1733427039589844, 10**18)
chan = {1: 1074, 2: 58, 3: 86, 4: 104, 5: 48, 6: 24, 7: 8, 8: 8}
sel = {1: 2412, 2: 418, 3: 432, 4: 142, 5: 82, 8: 8}
def add(b, e):
    o = dict(b)
    for r, k in e.items(): o[r] = o.get(r, 0) + k
    return o
mom = lambda loc, a, cen=56: CEIL.relax(loc, 73, 832, 9, cen, a)
g = Q(1, 10**12)
def maxint(pred, hi):
    lo = 0
    while hi - lo > 1:
        m = (lo+hi)//2
        try: ok = pred(m)
        except ValueError: ok = False
        lo, hi = (m, hi) if ok else (lo, m)
    return lo
D = dict(a_sel=float(asel), a_bit=float(abit),
         selected_ceiling=float(CEIL.ceiling(sel, 73, 896, 9, 64, grid=g)[1]),
         p104_channel_free_correction_ceiling=float(CEIL.ceiling(chan, 73, 832, 9, 56, grid=g)[1]),
         p104_minimal_word_ceiling=float(CEIL.ceiling({8: 208, 7: 8}, 73, 832, 9, 56, grid=g)[1]),
         idle_mass_per_call=str(CEIL.model(chan, 73, 832, 9, 56)['idle_mass']))
for nm, a in (('a0', Q(0)), ('a_sel', asel), ('a_bit', abit)):
    lo, hi = Q(0), Q(100)
    while hi - lo > Q(1, 10**6):
        mid = (lo+hi)/2
        lo, hi = (mid, hi) if mom(add(chan, {1: mid}), a, 56+mid)[1] < 1 else (lo, mid)
    D['nonproductive_rank1_units_max_' + nm] = round(float(lo), 4)
    D['productive_width9_helpers_max_' + nm] = maxint(lambda k: mom(add(chan, {1: 2*k, 7: k}), a)[1] < 1, 483)
D['moment_at_a_sel_with_96_nonproductive'] = float(mom(add(chan, {1: 96}), asel, 56+96)[0])
D['moment_at_0_with_96_nonproductive'] = str(Q(9*(2368+96) + 72*832 - 9*2312, 73*832))
D['ceiling_with_208_relay_helpers_117'] = float(CEIL.ceiling(add(chan, {1: 416, 7: 208}), 73, 832, 9, 56, grid=g)[1])
D['moment_at_a_sel_with_4480_noRelay_helpers'] = 'infeasible: 4480*9=40320 > idle mass ' + D['idle_mass_per_call']
R['D_budget'] = D
# ---------------- E. orbits (optional, slow-ish) ----------------
if '--skip-orbits' not in sys.argv:
    L = sorted({max(v, tuple(-x for x in v)) for v in rr}); N = len(L); idx = {l: i for i, l in enumerate(L)}
    nb = [[j for j in range(N) if j != i and ip(L[i], L[j]) == 0] for i in range(N)]
    adjs = [set(x) for x in nb]
    init = [1 if L[i] in removed else 0 for i in range(N)]
    def refine(col):
        while True:
            sig = [(col[i], tuple(sorted(collections.Counter(col[j] for j in nb[i]).items()))) for i in range(N)]
            keys = sorted(set(sig)); m = {s: k for k, s in enumerate(keys)}; new = [m[s] for s in sig]
            if len(keys) == len(set(col)): return new, keys
            col = new
    def extend(cA, cB):
        a, ka = refine(cA); b, kb = refine(cB)
        if ka != kb: return None
        if len(set(a)) == N:
            pos = {c: i for i, c in enumerate(b)}; perm = [pos[a[i]] for i in range(N)]
            return perm if all(perm[j] in adjs[perm[i]] for i in range(N) for j in nb[i]) else None
        cnt = collections.Counter(a); c = min((k for k in cnt if cnt[k] > 1), key=lambda k: cnt[k])
        u = a.index(c); top = max(a)+1; a2 = list(a); a2[u] = top
        for v in [i for i in range(N) if b[i] == c]:
            b2 = list(b); b2[v] = top; r_ = extend(a2, b2)
            if r_: return r_
        return None
    pidx = [idx[p] for p in ports]; parent = {i: i for i in pidx}
    def find(x):
        while parent[x] != x: x = parent[x]
        return x
    base, _ = refine(init); changed = True
    while changed:
        changed = False; reps = sorted({find(i) for i in pidx})
        for i in range(len(reps)):
            for j in range(i+1, len(reps)):
                x, y = reps[i], reps[j]
                if find(x) == find(y): continue
                a = list(base); b = list(base); top = max(base)+1; a[x] = top; b[y] = top
                p = extend(a, b)
                if p:
                    changed = True
                    for k in pidx:
                        s_, t_ = find(k), find(p[k])
                        if s_ != t_: parent[s_] = t_
    orb = collections.defaultdict(list)
    for i in pidx: orb[find(i)].append(ports.index(L[i]))
    orbs = sorted(orb.values(), key=len)
    R['E_orbits'] = dict(sizes=[len(o) for o in orbs],
        gamma_degree_matrix=[[sum(1 for i in a for j in b if S[i][j])//len(a) for b in orbs] for a in orbs],
        note='orbits of the stabiliser of F0uC in Aut(E8 line-orthogonality graph)=W(E8)/{+-1}; separations proven by exhaustive individualisation-refinement')
json.dump(R, open(os.path.join(HERE, 'RESULTS.json'), 'w'), indent=1, default=str)
print(json.dumps(R, indent=1, default=str)); print('PASS_EXPAZ_P104_CORRECTION_SCREEN')

"""Enumerate single-gate deferral/donor-substitution merges of every 0->1->2 rank-1 helper stop
in PR1122's word. A move: delete A gate g_i = (F1, h, [s]) and insert (F2, h, [r]) before global
time tau in (i, T2] where T2 = h's next touch, F2 = h's next frame, val_r(tau) == val_s(i)."""
import gzip, json, sys, importlib.util, bisect
from collections import Counter, defaultdict
ROOT = sys.argv[1]
W = json.loads(gzip.decompress(open(ROOT+'/research/r439-abc-semantic-bridge/inputs/candidate.json.gz','rb').read()))
spec = importlib.util.spec_from_file_location('xc', ROOT+'/research/e8-bit-p111-tuned/code/xcheck_p.py')
xc = importlib.util.module_from_spec(spec); spec.loader.exec_module(xc)
O = xc.Word(W)
v, R = W['v'], W['R']; n = 2*v+R
SC = len(W['D'])+len(W['A'])  # scatter time slot index (scatter happens before gate SC)
gates = [(ph,i,g) for ph in 'DABC' for i,g in enumerate(W[ph]) if ph!='C' or True]
# order: D, A, [scatter], B, C
gates = [('D',i,g) for i,g in enumerate(W['D'])]+[('A',i,g) for i,g in enumerate(W['A'])]+[('B',i,g) for i,g in enumerate(W['B'])]+[('C',i,g) for i,g in enumerate(W['C'])]
T = len(gates)
val = [1<<r for r in range(n)]
cur = list(W['start'])
touches = defaultdict(list)   # r -> [(t, frame)]
vhist = defaultdict(list)     # r -> [(t_from, value)]  value valid from before gate t_from
for r in range(n): vhist[r].append((0, val[r]))
for t,(ph,i,(f,d,src)) in enumerate(gates):
    if t == SC:
        copies = [val[q] for q in W['ret']]
        for tt,row in enumerate(W['scatter']):
            for k in row: val[v+tt] ^= copies[k]
        for tt in range(v): vhist[v+tt].append((t, val[v+tt]))
    for r in [d]+src: touches[r].append((t,f))
    acc = 0
    for s in src: acc ^= val[s]
    val[d] ^= acc
    vhist[d].append((t+1, val[d]))
def value_at(r, t):  # value before gate t
    h = vhist[r]; k = bisect.bisect_right([a for a,_ in h], t)-1
    return h[k][1]
def frame_before(r, t):  # frame register r holds just before gate t, and next frame it must reach after
    tl = touches[r]; k = bisect.bisect_left([a for a,_ in tl], t)
    curf = tl[k-1][1] if k>0 else W['start'][r]
    nxt = tl[k][1] if k < len(tl) else W['final'][r]
    return curf, nxt, (tl[k][0] if k < len(tl) else T)
def le(a,b):  # a subseteq b legal climb-or-stay
    return a==b or (O.nested(a,b) and O.fdim[a] < O.fdim[b])
# interval index of values
byval = defaultdict(list)
for r in range(n):
    h = vhist[r]
    for k,(a,x) in enumerate(h):
        b = h[k+1][0] if k+1 < len(h) else T+1
        byval[x].append((r,a,b))
# stops
paths = {}
stops = []
for r in range(2*v, n):
    tl = touches[r]
    fr = [W['start'][r]]; first_t = []
    for t,f in tl:
        if f != fr[-1]: fr.append(f); first_t.append(t)
    if len(fr) >= 3 and O.fdim[fr[0]]==0 and O.fdim[fr[1]]==1 and O.fdim[fr[2]]==2:
        at1 = [t for t,f in tl if f==fr[1]]
        stops.append((r, fr[1], fr[2], at1, first_t[1]))
print('stops', len(stops))
res = Counter(); legal = []; witness = {}
for h,F1,F2,at1,T2 in stops:
    assert len(at1)==1
    t0 = at1[0]; ph,i,(f,d,src) = gates[t0]
    assert ph=='A' and d==h and len(src)==1, gates[t0]
    s = src[0]; X = value_at(s, t0)
    reasons = Counter(); found = []
    for (r,a,b) in byval[X]:
        if r == h or v <= r < 2*v: continue
        lo, hi = max(a, t0+1 if r!=s else t0), min(b, T2)  # tau in [lo, hi]; tau=t0 for r=s means reframe in place
        # for r==s: value valid at t0 too
        if r == s: lo = max(a, t0)
        if lo > hi: continue
        tau = lo
        while tau <= hi:
            curf, nxt, nt = frame_before(r, tau)
            # careful: if r==s and tau<=t0 the original gate is still in the list -> treat touches of s at t0 as removed
            if r == s and nt == t0:
                tl = touches[r]; k = [a for a,_ in tl].index(t0)
                nxt = tl[k+1][1] if k+1 < len(tl) else W['final'][r]; nt2 = tl[k+1][0] if k+1<len(tl) else T
            else: nt2 = nt
            ok1, ok2 = le(curf, F2), le(F2, nxt)
            if ok1 and ok2: found.append((r, tau, curf, nxt))
            else: reasons[(('cur' if not ok1 else '')+('next' if not ok2 else '')) or 'x'] += 1
            tau = nt2 + 1 if nt2 >= tau else hi+1
            if nt2 >= hi: break
    if found: legal.append((h, t0, F2, found)); res['has_candidate'] += 1
    else: res['none'] += 1
    witness[h] = dict(stop_gate=['A', i, [f,d,src]], F1=W['frames'][F1], F2=W['frames'][F2], F2id=F2, T2=T2,
                      donors_with_value=sorted({r for (r,a,b) in byval[X] if a<=T2 and b>t0 and r!=h}),
                      reasons=dict(reasons), source_path_frames=[W['frames'][fr] for _,fr in touches[s] if True][:0])
print(res)
for L in legal[:20]: print(L[:3], L[3][:6])
json.dump({'stops':len(stops),'legal':[[a,b,c,d] for a,b,c,d in legal],'witness':witness}, open(sys.argv[2],'w'))

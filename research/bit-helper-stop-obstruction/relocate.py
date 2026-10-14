"""Generalized single-gate relocation of each rank-1 stop gate g=(F1,h,[s]) of PR1122:
delete g, insert (F,h,[r]) before global time tau, t0 < tau <= first READ of h after t0,
val_r(tau)==val_s(t0), r not a target, F any closed frame (catalogue or closure) with
cur_h(tau) <= F <= next_h(tau) and cur_r(tau) <= F <= next_r(tau) (<= : equal or nested strict climb).
Recompute h's and r's paths, profile delta and moment delta  sum n_d d (64/d)^c."""
import gzip, json, sys, importlib.util, bisect
from collections import Counter, defaultdict
from fractions import Fraction
ROOT = sys.argv[1]
W = json.loads(gzip.decompress(open(ROOT+'/research/r439-abc-semantic-bridge/inputs/candidate.json.gz','rb').read()))
spec = importlib.util.spec_from_file_location('xc', ROOT+'/research/e8-bit-p111-tuned/code/xcheck_p.py')
xc = importlib.util.module_from_spec(spec); spec.loader.exec_module(xc)
O = xc.Word(W)
v, R = W['v'], W['R']; n = 2*v+R
FR = W['frames']; fset = [frozenset(f) for f in FR]
SC = len(W['D'])+len(W['A'])
gates = [('D',i,g) for i,g in enumerate(W['D'])]+[('A',i,g) for i,g in enumerate(W['A'])]+[('B',i,g) for i,g in enumerate(W['B'])]+[('C',i,g) for i,g in enumerate(W['C'])]
T = len(gates)
val = [1<<r for r in range(n)]
touches = defaultdict(list); vhist = defaultdict(list)
for r in range(n): vhist[r].append((0, val[r]))
for t,(ph,i,(f,d,src)) in enumerate(gates):
    if t == SC:
        copies = [val[q] for q in W['ret']]
        for tt,row in enumerate(W['scatter']):
            for k in row: val[v+tt] ^= copies[k]
        for tt in range(v): vhist[v+tt].append((t, val[v+tt]))
    for r in [d]+src: touches[r].append((t,f,'d' if r==d else 's'))
    acc = 0
    for s in src: acc ^= val[s]
    val[d] ^= acc
    vhist[d].append((t+1, val[d]))
def le(a,b): return a==b or (O.nested(a,b) and O.fdim[a] < O.fdim[b])
def dim(f): return O.fdim[f]
def path_of(r, tl):
    p=[W['start'][r]]
    for t,f,_ in tl:
        if f!=p[-1]: p.append(f)
    if W['final'][r]!=p[-1]: p.append(W['final'][r])
    return p
def prof(p): return Counter(dim(b)-dim(a) for a,b in zip(p,p[1:]))
c = Fraction(1733427039589844234749, 10**24)
def cost(cn):  # float display of sum n_d d (64/d)^c
    return sum(k*d*(64/d)**float(c) for d,k in cn.items())
MASK = ((1<<(2*v))-1) if '--dirty' in sys.argv else -1
byval = defaultdict(list)
for r in range(n):
    h = vhist[r]
    for k,(a,x) in enumerate(h):
        b = h[k+1][0] if k+1 < len(h) else T+1
        byval[x & MASK].append((r,a,b))
# all frames in catalogue as candidates (636), plus none extra: closures of unions checked separately
stops=[]
for h in range(2*v,n):
    p = path_of(h, touches[h])
    if len(p)>=3 and dim(p[0])==0 and dim(p[1])==1 and dim(p[2])==2: stops.append(h)
print('stops',len(stops))
best=[]; stats=Counter(); per={}
for h in stops:
    tl = touches[h]; p1 = path_of(h, tl)[1]
    k1 = [x[1] for x in tl].index(p1); assert sum(1 for x in tl if x[1]==p1)==1
    t0,f0,_ = tl[k1]; ph,i,(f,d,src) = gates[t0]
    assert d==h and len(src)==1 and dim(f)==1
    s = src[0]; X = vhist[s][bisect.bisect_right([a for a,_ in vhist[s]], t0)-1][1]
    # first read of h after t0
    tread = next((t for t,ff,role in tl if t>t0 and role=='s'), T)
    tl_h = [x for x in tl if x[0]!=t0]
    base_h = prof(path_of(h, tl)); ncand=0; nlegal=0; bestd=None
    for (r,a,b) in byval[X & MASK]:
        if r==h or v<=r<2*v: continue
        lo, hi = max(a, t0), min(b-1, tread)   # tau in [lo,hi]; value valid before gate tau for tau in [a,b)
        if lo>hi: continue
        tlr = [x for x in touches[r] if x[0]!=t0] if r==s else touches[r]
        times_r = [x[0] for x in tlr]; times_h=[x[0] for x in tl_h]
        # candidate taus: breakpoints of either register's touch sequence
        taus = sorted({lo,hi} | {t for t in times_r+times_h if lo<=t<=hi} | {t+1 for t in times_r+times_h if lo<=t+1<=hi})
        for tau in taus:
            ncand += 1
            def cn(tlx, reg):
                k = bisect.bisect_left([x[0] for x in tlx], tau)
                cur = tlx[k-1][1] if k>0 else W['start'][reg]
                nxt = tlx[k][1] if k<len(tlx) else W['final'][reg]
                return k,cur,nxt
            kh,ch,nh = cn(tl_h,h); kr,cr,nr = cn(tlr,r)
            for F in range(len(FR)):
                if not (le(ch,F) and le(F,nh) and le(cr,F) and le(F,nr)): continue
                nlegal += 1
                newh = tl_h[:kh]+[(tau,F,'d')]+tl_h[kh:]
                if r==s:
                    newr = tlr[:kr]+[(tau,F,'s')]+tlr[kr:]
                    oldr = touches[r]
                else:
                    newr = tlr[:kr]+[(tau,F,'s')]+tlr[kr:]; oldr = touches[r]
                delta = prof(path_of(h,newh)); delta.subtract(base_h)
                dr = prof(path_of(r,newr)); dr.subtract(prof(path_of(r,oldr))); delta.update(dr)
                if r!=s:
                    ds = prof(path_of(s,[x for x in touches[s] if x[0]!=t0])); ds.subtract(prof(path_of(s,touches[s]))); delta.update(ds)
                dd = {k:x for k,x in delta.items() if x}
                dc = cost(Counter({k:x for k,x in dd.items()}))
                if bestd is None or dc < bestd[0]: bestd=(dc, dd, r, tau, F)
                if dc < -1e-12: best.append((h,t0,r,tau,F,dd,dc))
    per[h]=dict(donors=sorted({r for (r,a,b) in byval[X & MASK] if r!=h and not v<=r<2*v and max(a,t0)<=min(b-1,tread)}),candidates=ncand, legal_frame_choices=nlegal, best=bestd and [bestd[0], bestd[1], bestd[2], bestd[3], FR[bestd[4]]])
    stats['legal' if nlegal else 'nolegal'] += 1
print(stats, 'improving', len(best))
print(Counter(json.dumps(x['best'][1] if x['best'] else None, sort_keys=True) for x in per.values()))
for b in best[:10]: print(b)
json.dump({'per':per,'improving':[list(map(str,b)) for b in best]}, open(sys.argv[2],'w'))

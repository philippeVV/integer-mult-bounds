import gzip, json, sys, importlib.util
from collections import Counter, defaultdict
ROOT=sys.argv[1]
W = json.loads(gzip.decompress(open(ROOT+'/research/r439-abc-semantic-bridge/inputs/candidate.json.gz','rb').read()))
spec = importlib.util.spec_from_file_location('xc', ROOT+'/research/e8-bit-p111-tuned/code/xcheck_p.py')
xc = importlib.util.module_from_spec(spec); spec.loader.exec_module(xc); O = xc.Word(W)
v,R=W['v'],W['R']; n=2*v+R; FR=W['frames']
gates = [('D',i,g) for i,g in enumerate(W['D'])]+[('A',i,g) for i,g in enumerate(W['A'])]+[('B',i,g) for i,g in enumerate(W['B'])]+[('C',i,g) for i,g in enumerate(W['C'])]
touch=defaultdict(list)
for t,(ph,i,(f,d,src)) in enumerate(gates):
    for r in [d]+src: touch[r].append((t,f,'d' if r==d else 's'))
def segs(r):
    out=[]; cur=W['start'][r]; first=None
    for t,f,role in touch[r]:
        if f!=cur or first is None:
            if f!=cur: out.append((cur, first, last)); cur=f; first=t
            elif first is None: first=t
        last=t
    out.append((cur,first,last)); return out
stops=[]
for h in range(2*v,n):
    p=[W['start'][h]]
    for t,f,_ in touch[h]:
        if f!=p[-1]: p.append(f)
    if len(p)>=3 and O.fdim[p[1]]==1 and O.fdim[p[2]]==2:
        t0=[t for t,f,_ in touch[h] if f==p[1]][0]; T2=[t for t,f,_ in touch[h] if f==p[2]][0]
        s=gates[t0][2][2][0]
        tread=next((t for t,f,role in touch[h] if t>t0 and role=='s'),None)
        stops.append((h,s,p[2],t0,T2,tread))
bys=defaultdict(list)
for x in stops: bys[x[1]].append(x)
print('sources with stops', len(bys), Counter(len(x) for x in bys.values()))
cls=Counter(); ex={}
for s,L in sorted(bys.items()):
    sg=[(FR[f],a,b) for f,a,b in segs(s)]
    # x_s's frame right after its rank-1 phase
    s_frames=[f for f,_,_ in segs(s)]
    for (h,_,F2,t0,T2,tread) in L:
        contains=[f for f in s_frames if O.nested(F2,f)]
        first_sup=contains[0] if contains else None
        # time x_s first sits on a frame containing F2
        tsup=next((a for f,a,b in segs(s) if O.nested(F2,f)), None)
        key=('s_never_on_F2' if F2 not in s_frames else 's_visits_F2',
             'first_superframe_after_h_read' if tsup is not None and tread is not None and tsup>tread else 'first_superframe_before_or_no_read',
             'written_between' if any(role=='d' and t0<t<(tread or 10**9) for t,f,role in touch[s]) else 'unwritten')
        cls[key]+=1
        if key not in ex: ex[key]=dict(h=h,s=s,stop_gate=gates[t0],F2=FR[F2],T2=T2,first_read_of_h=tread,s_segments=sg[:6],tsup=tsup,first_superframe_dim=O.fdim[first_sup] if first_sup is not None else None)
print(cls)
for k,e in ex.items(): print(k, e)
print('--- late rank-1 touch witness')
wc=Counter(); rows=[]
for (h,s,F2,t0,T2,tread) in stops:
    r1=[(t,role) for t,f,role in touch[s] if O.fdim[f]==1]
    last=max(t for t,_ in r1)
    g=gates[last]
    kind=(g[0], 'tgt' if v<=g[2][1]<2*v else ('src' if g[2][1]<v else 'hlp'), 'x_s_is_dest' if g[2][1]==s else 'x_s_donor')
    ok = tread is not None and last > tread
    wc[(ok,)+kind]+=1
    rows.append(dict(h=h,s=s,F2=FR[F2],stop_gate=['A',gates[t0][1],gates[t0][2]],h_first_read=tread,h_first_read_gate=[gates[tread][0],gates[tread][1]],s_last_rank1_touch=last,s_last_rank1_gate=[g[0],g[1],g[2]]))
print(wc)
json.dump(rows, open(sys.argv[2],'w'), indent=0)
print('--- refined')
rc=Counter()
for (h,s,F2,t0,T2,tread) in stops:
    r1=[t for t,f,role in touch[s] if O.fdim[f]==1]; last=max(r1)
    nxt=next((t,f) for t,f,role in touch[s] if t>last)
    tG,G=nxt
    if last>tread: k='A: x_s still has rank-1 touch at {s} after h first read (B target write)'
    elif not O.nested(F2,G): k='B: x_s next frame G does not contain F2 (dim %d)'%O.fdim[G]
    else: k='C: G contains F2 (dim %d); inserting F2 stop on x_s splits 1->%d'%(O.fdim[G],O.fdim[G])
    rc[k]+=1
for k,c in sorted(rc.items()): print(c,k)

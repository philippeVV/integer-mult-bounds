"""Hostee-move / augmenting census on a BORROW-WORD.
lives: B-warp hostees (record k: register X, life = X's events after the warp) and free followers (helpers entering in B).
registers: own life = events before the first B-warp on it. edge(life->reg Y): Y not restored, own-life last < life first,
dim lo(Y own last gate) <= dim hi(life first gate)  [refit-necessary, lattice bound], and (for followers) frame rule as pairs.py.
Prints for each free follower the augmenting paths (depth<=3) that host it while keeping every hostee hosted."""
import sys,json
sys.path.insert(0,sys.argv[2])
from collections import defaultdict,Counter
from bounds import Bounds,regs
b=Bounds(sys.argv[1]);Wd=b.Wd;w=Wd.w;v=Wd.v;ev=b.ev
restored={int(q) for q in w.get('helper_final',{})}
touch=defaultdict(list)
for t,(sec,e) in enumerate(ev):
    for q in ((e[1],) if e[0]=='warp' else regs(e)):touch[q].append((t,sec,e[0]))
def own_and_hosted(r):
    T=touch[r];bw=[i for i,(t,sec,k) in enumerate(T) if k=='warp' and sec=='Bwarp']
    if not bw:return T,[]
    return T[:bw[0]],T[bw[0]:]
lives={};regsinfo={}
for r in Wd.live:
    if v<=r<2*v:continue
    own,host=own_and_hosted(r)
    if host:
        g=[t for t,sec,k in host if k!='warp']
        lives[('H',r)]=dict(first=host[0][0],firstg=g[0] if g else None,hi=len(b.hi[g[0]]) if g else 9,
                            restore=any(sec=='restore' for _,sec,_ in host),cur=r,n=len(g))
    if r in restored:continue
    if any(sec=='restore' for _,sec,_ in own):continue
    lastf=[x for x in b.byreg.get(r,[]) if x[0]<=(own[-1][0] if own else -1)]
    if lastf and lastf[-1][1]=='gate':lo=len(b.lo[lastf[-1][0]])
    elif lastf:lo=len(lastf[-1][3])
    else:lo=len(Wd.start[r])
    regsinfo[r]=dict(last=own[-1][0] if own else -1,lo=lo,hosted=bool(host))
lenders_now={rec[1] for rec in w['warps']}
for h in Wd.helpers:
    T=touch.get(h,[])
    if h in restored or not T or any(sec not in('B','restore','Bwarp') for _,sec,_ in T) or h in lenders_now:continue
    g=[t for t,sec,k in T if k!='warp']
    lives[('F',h)]=dict(first=g[0],firstg=g[0],hi=len(b.hi[g[0]]),cur=None,end=T[-1][0],restore=any(sec=='restore' for _,sec,_ in T),n=len(g))
def ok(L,Y):
    l=lives[L];y=regsinfo.get(Y)
    return y is not None and y['last']<l['first'] and y['lo']<=l['hi']
occ={l['cur']:L for L,l in lives.items() if l['cur'] is not None}
F=[L for L in lives if L[0]=='F']
print('followers',[(L,lives[L]['first'],lives[L]['hi']) for L in F if lives[L]['hi']>=4])
print('hostees',len([L for L in lives if L[0]=='H']),'reg lo hist',Counter((r['lo'],r['hosted']) for r in regsinfo.values()))
# direct edges summary for hostees: alternative registers
for L,l in sorted(lives.items()):
    if L[0]!='H':continue
    alt=[Y for Y in regsinfo if Y!=l['cur'] and ok(L,Y)]
    freealt=[Y for Y in alt if Y not in occ]
    print('hostee',L,'first',l['first'],'hi',l['hi'],'restore',l['restore'],'curlo',regsinfo.get(l['cur'],{}).get('lo'),'free alts',freealt[:12],'n',len(freealt),'occ alts',len(alt)-len(freealt))
# augmenting paths for followers with hi>=4
def aug(L,used,depth):
    for Y in regsinfo:
        if Y in used or not ok(L,Y):continue
        if Y not in occ:yield [(L,Y)]
        elif depth>0:
            for p in aug(occ[Y],used|{Y},depth-1):yield [(L,Y)]+p
for L in F:
    if lives[L]['hi']<4:continue
    P=list(aug(L,set(),3));print('follower',L,'paths',len(P));[print('  ',p) for p in P[:40]]
print('ends: followers',[(L,lives[L]['first'],lives[L]['end'],lives[L]['n']) for L in F if lives[L]['hi']>=4])
he=sorted((touch[L[1]][-1][0],L) for L in lives if L[0]=='H');print('hostee life ends min',he[:5],'max',he[-1],'nevents',len(ev))
# any register hosting a life that ends before another life's first? (sequential double hosting)
fol_end={L:lives[L]['end'] for L in F}
fr=sorted(lives[L]['first'] for L in lives if L[0]=='H')
print('earliest hostee first',fr[:6])

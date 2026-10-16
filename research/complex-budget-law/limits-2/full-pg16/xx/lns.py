"""CP-SAT large-neighbourhood search on exact F: pick a vertex set X (front+back blocks across the seam, a contiguous window,
or a random scatter); CP-SAT re-permutes X among X's slots and re-optimises the orientation of ALL edges (in-degrees 8/9),
minimising #back arcs exactly.  usage: lns.py orderfile rounds size secs seed out"""
import sys, random, numpy as np
from ortools.sat.python import cp_model
from common import *
order=[int(x) for x in open(sys.argv[1]).read().split()][:V]
rounds,size,secs,seed,out=int(sys.argv[2]),int(sys.argv[3]),float(sys.argv[4]),int(sys.argv[5]),sys.argv[6]
rng=random.Random(seed)
import os
KINDS=os.environ.get('KINDS','seam,window,scatter,seam').split(',')
cur=F(order); print('start F',cur,flush=True)
def solve(order,slots):
    pos=np.empty(V,int); pos[np.asarray(order)]=np.arange(V)
    X=[order[s] for s in slots]; inX={v:t for t,v in enumerate(X)}; n=len(X)
    m=cp_model.CpModel()
    p={v:m.NewIntVar(0,n-1,f'p{v}') for v in X}; m.AddAllDifferent(list(p.values()))
    slots_arr=np.array(slots)
    o=[m.NewBoolVar(f'o{e}') for e in range(len(ED))]     # o=1: arc ED[e][0]->ED[e][1] (point->line)
    indeg=[[] for _ in range(V)]
    back=[]
    cache={}
    for e,(a,b) in enumerate(ED):
        indeg[b].append(o[e]); indeg[a].append(o[e].Not())
        if a not in inX and b not in inX:
            bef=pos[a]<pos[b]          # a before b
            back.append(o[e].Not() if bef else o[e]); continue
        # before_ab boolean
        if a in inX and b in inX:
            key=(a,b); bb=m.NewBoolVar(''); m.Add(p[a]<p[b]).OnlyEnforceIf(bb); m.Add(p[a]>p[b]).OnlyEnforceIf(bb.Not())
        else:
            x,u=(a,b) if a in inX else (b,a)
            c=int(np.searchsorted(slots_arr,pos[u]))      # x before u iff p_x < c
            if c==0: xb=None; val=False
            elif c==n: xb=None; val=True
            else:
                xb=m.NewBoolVar(''); m.Add(p[x]<c).OnlyEnforceIf(xb); m.Add(p[x]>=c).OnlyEnforceIf(xb.Not())
            if xb is None:
                bef_xu=val; bef_ab=bef_xu if x==a else (not bef_xu)
                back.append(o[e].Not() if bef_ab else o[e]); continue
            bb=xb if x==a else xb.Not()
        # back iff o != bb  : arc a->b with b before a, or arc b->a with a before b
        bk=m.NewBoolVar(''); m.AddBoolXOr([bk,o[e],bb.Not()]) if False else None
        # bk == (o XOR bb)
        m.Add(bk>=o[e]-bb); m.Add(bk>=bb-o[e]); m.Add(bk<=o[e]+bb); m.Add(bk<=2-o[e]-bb)
        back.append(bk)
    for v in range(V): m.Add(sum(indeg[v])==int(K[v]))
    m.Minimize(sum(back))
    # hint: current
    for v in X: m.AddHint(p[v],inX[v])
    _,t_,h_,_=F(order,True); dirmap={(int(a),int(b)) for a,b in zip(t_,h_)}
    for e,(a,b) in enumerate(ED): m.AddHint(o[e], 1 if (int(a),int(b)) in dirmap else 0)
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=secs; s.parameters.num_workers=2
    st=s.Solve(m)
    if st not in (cp_model.OPTIMAL,cp_model.FEASIBLE): return None,None
    newX=sorted(X,key=lambda v:s.Value(p[v])); new=list(order)
    for t,sl in enumerate(slots): new[sl]=newX[t]
    global LAST; LAST=(s.StatusName(st),s.BestObjectiveBound())
    return int(s.ObjectiveValue()),new
best=cur
for r in range(rounds):
    kind=rng.choice(KINDS)
    if kind=='seam':
        a=rng.randint(size//4,3*size//4); slots=list(range(a))+list(range(V-(size-a),V))
    elif kind=='window':
        st=rng.randrange(V-size); slots=list(range(st,st+size))
    else:
        slots=sorted(rng.sample(range(V),size))
    val,new=solve(order,slots)
    if val is None: print(r,kind,'no sol',flush=True); continue
    Fn=F(new)
    if Fn<=cur:
        order=new
        if Fn<cur: print(f'round {r} {kind}: F {cur} -> {Fn}',flush=True)
        cur=Fn
        if cur<best:
            best=cur; open(out,'w').write('\n'.join(map(str,order))+'\n')
    if True: print(f'round {r} {kind} cp {val} {LAST} F {Fn} cur {cur} best {best}',flush=True)
print('final best',best)

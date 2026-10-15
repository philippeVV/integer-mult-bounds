import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""SA v2 on the hub relaxation, flags stored as subspace lists; 'reroute' moves snap a register's chain
through levels already used by other registers (sharing is what moves content).
Registers: 112 targets (L1=<p_sigma(t)>, L8=p_t^perp), 112 sources (L1=<p_s>, L2=<p_s,u_s>, u_s in K_s,
L8 = p_tau(s)^perp); sigma,tau bijections (perfect matchings); --lost L sources get no top (centre slots);
--k helpers: chain from a start pin of dim (8-width) to a top, cost f(width)."""
import json, random, math, sys, argparse, time
sys.path.insert(0,'../expBV'); from wordmodel import basis
ap=argparse.ArgumentParser(); ap.add_argument('--iters',type=int,default=100000); ap.add_argument('--seed',type=int,default=1)
ap.add_argument('--lost',type=int,default=0); ap.add_argument('--k',type=int,default=0); ap.add_argument('--hw',type=int,default=2)
ap.add_argument('--T0',type=float,default=4.0); ap.add_argument('--out',default=''); ap.add_argument('--ufree',action='store_true')
a=ap.parse_args(); random.seed(a.seed)
P=json.load(open('instA.json'))['ports']; cen=json.load(open('cen.json'))['cen']; n=112
nor=[1,7,25,31,121,127,97,103]; dot=lambda x,y:(x&y).bit_count()&1
PERP=[basis([x for x in range(1,512) if dot(x,P[t])==0]) for t in range(n)]
K={s:[x for x in range(1,512) if x!=P[s] and (a.ufree or all(dot(x,nor[q])==0 for q in cen[s]))] for s in range(n)}
def sp(B):
    s=[0]
    for v in B: s=s+[x^v for x in s]
    return s
def contained(A,B): return basis(list(A)+list(B))==B
def between(A,B):
    """all subspaces C with A < C < B, dim C = dim A+1"""
    out=set()
    for x in sp(B):
        if x and not contained((x,),A) : out.add(basis(list(A)+[x]))
    return list(out)
def randflag(lo_key, top, mid=None):
    """random full flag from lo_key up to top"""
    ch=[lo_key]
    while len(ch[-1])<len(top)-1:
        ch.append(random.choice(between(ch[-1],top)))
    return ch+[top]
# perfect matchings sigma (target->source) and tau (source->top) with constraints
def matching(adj, nL):
    m={}
    def tryk(u,seen):
        for v in random.sample(adj[u],len(adj[u])):
            if v in seen: continue
            seen.add(v)
            if v not in m or tryk(m[v],seen): m[v]=u; return True
        return False
    for u in range(nL): tryk(u,set())
    return {u:v for v,u in m.items()}
sig=matching({t:[s for s in range(n) if dot(P[s],P[t])==0] for t in range(n)},n)
U={s:random.choice([u for u in K[s] if any(dot(P[t],P[s])==0 and dot(P[t],u)==0 for t in range(n))]) for s in range(n)}
lost=set(random.sample(range(n),a.lost))
tau=matching({s:[t for t in range(n) if dot(P[t],P[s])==0 and dot(P[t],U[s])==0] if s not in lost else [] for s in range(n)},n)
print('matched sigma',len(sig),'tau',len(tau),flush=True)
regs=[]; FL=[]
for t in range(n): regs.append(('T',t,1)); FL.append([(P[sig[t]],)]+randflag((P[sig[t]],),PERP[t])[1:] if True else None)
for s in range(n):
    if s in tau:
        G=basis([P[s],U[s]]); regs.append(('S',s,2)); FL.append([(P[s],)]+randflag(G,PERP[tau[s]]))
    else:
        regs.append(('X',s,9)); FL.append(None)   # lost: no top, ignore post-cut (favourable: its content still via island)
for i in range(a.k):
    t=random.randrange(n); lo=8-a.hw; ch=randflag(basis([random.choice(PERP[t])]),PERP[t]); regs.append(('H',i,lo)); FL.append(ch)
NR=len(regs)
for i in range(NR):
    if FL[i] is not None: assert len(FL[i])==8 and [len(x) for x in FL[i]]==list(range(1,9)), (regs[i],[len(x) for x in FL[i]])
base={(P[s],):1<<s for s in range(n)}
def score():
    cont=[0]*NR; top=None
    for d in range(8):
        g={}
        for i in range(NR):
            f=FL[i]
            if f is None: continue
            if regs[i][0]=='H' and d<regs[i][2]-1: continue
            kk=f[d]; g[kk]=g.get(kk,0)|cont[i]
        for kk in g: g[kk]|=base.get(kk,0)
        for i in range(NR):
            f=FL[i]
            if f is None or (regs[i][0]=='H' and d<regs[i][2]-1): continue
            cont[i]=g[f[d]]
    tot=0
    for t in range(n): tot+=bin(cont[t]).count('1')
    return tot, cont
def levels_at(d):
    s={}
    for i in range(NR):
        if FL[i] is not None: s[FL[i][d]]=s.get(FL[i][d],0)+1
    return s
cur,_=score(); best=cur; bestFL=None; t0=time.time()
for it in range(a.iters):
    T=a.T0*(1-it/a.iters)+0.05
    i=random.randrange(NR)
    if FL[i] is None: continue
    kind,_,lo=regs[i]; old=list(FL[i])
    lo=max(lo,1)  # fixed levels: indices < lo are fixed (T: L1 fixed -> lo=1; S: L1,L2 fixed -> lo=2)
    if kind=='H': lo=0 if False else lo-1
    if lo>=7: continue
    r=random.random()
    if r<0.5:
        d=random.randrange(lo,7)
        FL[i][d]=random.choice(between(FL[i][d-1],FL[i][d+1])) if d>0 else FL[i][d]
    else:
        a0=random.randrange(max(lo-1,0),7); b0=random.randrange(a0+2,9) if a0+2<=8 else 8
        b0=min(b0,7)
        if b0<=a0+1: continue
        ch=FL[i]
        for d in range(a0+1,b0):
            cands=[L for L in levels_at(d) if L!=ch[d] and contained(ch[d-1],L) and contained(L,ch[b0])]
            if cands and random.random()<0.9: ch[d]=random.choice(cands)
            else: ch[d]=random.choice(between(ch[d-1],ch[b0])) if len(ch[b0])-len(ch[d-1])>=2 else ch[d]
            # keep nesting valid for d+1..b0
        ok=all(contained(ch[d],ch[d+1]) for d in range(7))
        if not ok: FL[i]=old; continue
    new,_=score()
    if new>=cur or random.random()<math.exp((new-cur)/T):
        cur=new
        if cur>best: best=cur; bestFL=[list(f) if f else None for f in FL]
    else: FL[i]=old
    if it%5000==0: print(it,cur,best,round(time.time()-t0),flush=True)
print('final best',best,'of 5824')
if a.out: json.dump(dict(best=best,regs=regs,FL=bestFL,args=vars(a)),open(a.out,'w'))

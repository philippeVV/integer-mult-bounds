import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""Simulated annealing on the hub (time-free) relaxation of the post-cut relay network.
Every register = full flag of port spans; constant hubs merge content at every shared level;
content x_s is born at level <p_s>. Targets: L1=<p_sigma(t)>, L8=p_t^perp. Sources: L1=<p_s>, L2=<p_s,u_s>
(u_s in K_s unless --ufree), L8 = p_tau(s)^perp (a target top) unless the source is 'lost' (centre slot).
Extra helpers (--k): flag from a pinned start level to a target top; cost = width = 8 - dim(start pin).
Score = #(s,t) pairs delivered (max 5824)."""
import json, random, math, sys, argparse
sys.path.insert(0,'../expBV'); from wordmodel import basis
ap=argparse.ArgumentParser(); ap.add_argument('--iters',type=int,default=200000); ap.add_argument('--seed',type=int,default=1)
ap.add_argument('--ufree',action='store_true'); ap.add_argument('--lost',type=int,default=0); ap.add_argument('--k',type=int,default=0)
ap.add_argument('--hstart',type=int,default=7); ap.add_argument('--T0',type=float,default=3.0); ap.add_argument('--out',default='')
a=ap.parse_args(); random.seed(a.seed)
P=json.load(open('instA.json'))['ports']; cen=json.load(open('cen.json'))['cen']; n=112
nor=[1,7,25,31,121,127,97,103]; dot=lambda x,y:(x&y).bit_count()&1
perp={t:[x for x in range(1,512) if dot(x,P[t])==0] for t in range(n)}
K={s:[x for x in range(1,512) if all(dot(x,nor[q])==0 for q in cen[s])] for s in range(n)}
orth={t:[s for s in range(n) if dot(P[s],P[t])==0] for t in range(n)}
def rank(vs): return len(basis(vs))
def complete(pref, top, old):
    """extend independent prefix to a basis of p_top^perp, prefer vectors from old, then random"""
    out=list(pref)
    for x in old+random.sample(perp[top],len(perp[top])):
        if len(out)==8: break
        if dot(x,P[top])==0 and rank(out+[x])==len(out)+1: out.append(x)
    assert len(out)==8 and all(dot(x,P[top])==0 for x in out); return out
# registers: ('T',t) ('S',s) ('H',i)
regs=[]; F=[]; info=[]
for t in range(n):
    s=random.choice(orth[t]); regs.append(('T',t)); info.append(dict(top=t)); F.append(complete([P[s]],t,[]))
lost=set(random.sample(range(n),a.lost)) if a.lost else set()
for s in range(n):
    regs.append(('S',s))
    us=[u for u in (range(1,512) if a.ufree else K[s]) if u!=P[s]]
    for _ in range(1000):
        u=random.choice(us); tops=[t for t in range(n) if dot(P[t],P[s])==0 and dot(P[t],u)==0]
        if tops: break
    t=random.choice(tops); info.append(dict(top=t,lost=s in lost)); F.append(complete([P[s],u],t,[]))
for i in range(a.k):
    t=random.randrange(n); regs.append(('H',i)); info.append(dict(top=t)); F.append(complete([],t,[]))
NR=len(regs)
def lo(i):
    k=regs[i][0]; return 2 if k=='T' else (3 if k=='S' else a.hstart+1)
def keys(i):
    f=F[i]; return [basis(f[:d]) for d in range(1,9)]
KEY=[keys(i) for i in range(NR)]
base={}
for s in range(n): base[(P[s],)]=base.get((P[s],),0)|(1<<s)
def score():
    cont=[0]*NR
    for d in range(1,9):
        g={}
        for i in range(NR):
            if regs[i][0]=='H' and d<=a.hstart: continue   # helper sits at its start pin until dim hstart
            if info[i].get('lost') and d==8: continue
            kk=KEY[i][d-1]; g[kk]=g.get(kk,0)|cont[i]
        for kk in g:
            if d==1: g[kk]|=base.get(kk,0)
        for i in range(NR):
            if regs[i][0]=='H' and d<a.hstart: continue
            if info[i].get('lost') and d==8: continue
            kk=KEY[i][d-1]
            cont[i]=g.get(kk,0)|base.get(kk,0) if d>1 or True else 0
            if d==1 and kk in base: cont[i]|=base[kk]
        if d==8: top=g
    return sum(bin(top.get(KEY[t][7],0)).count('1') for t in range(n)), top
cur,_=score(); best=cur; bestF=None
T0=a.T0
for it in range(a.iters):
    T=T0*(1-it/a.iters)+1e-3
    i=random.randrange(NR); kind=regs[i][0]; old=list(F[i]); oldinfo=dict(info[i])
    r=random.random()
    if r<0.8:
        d=random.randrange(lo(i),8)  # change level d (1-based index d means L_d, swap b_d with b_{d+1} or add)
        f=F[i]
        if random.random()<0.5: f[d-1],f[d]=f[d],f[d-1]
        else: f[d-1]^=f[d]
    elif kind=='T':
        s=random.choice(orth[regs[i][1]]); F[i]=complete([P[s]],info[i]['top'],old[1:])
    elif kind=='S':
        s=regs[i][1]; u=old[1]
        if random.random()<0.5:
            us=[x for x in (range(1,512) if a.ufree else K[s]) if x!=P[s]]; u=random.choice(us)
        tops=[t for t in range(n) if dot(P[t],P[s])==0 and dot(P[t],u)==0]
        if not tops: continue
        t=random.choice(tops); info[i]['top']=t; F[i]=complete([P[s],u],t,old[2:])
    else:
        t=random.randrange(n); info[i]['top']=t; F[i]=complete([],t,old)
    KEY[i]=keys(i); new,_=score()
    if new>=cur or random.random()<math.exp((new-cur)/T):
        cur=new
        if cur>best: best=cur; bestF=[list(f) for f in F]
    else:
        F[i]=old; info[i]=oldinfo; KEY[i]=keys(i)
    if it%20000==0: print(it,cur,best,flush=True)
print('final best',best,'of 5824')
if a.out: json.dump(dict(best=best,F=bestF,regs=regs,args=vars(a)),open(a.out,'w'))

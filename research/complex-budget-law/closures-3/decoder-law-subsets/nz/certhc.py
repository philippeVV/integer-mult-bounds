import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
# Lemma-B / L8 certificate ceiling on the frame-gauge P120 decoder:
# maximise rank_F2 (A Bh^T)^o over 0/1 labels A <= supp P, Bh <= supp Q (o = keep diag+orthogonal pairs)
import random, sys
from e8 import *
fr = frames()[0]
M = [ROOT[f] for f in fr]
terms = gauge_terms(M)                     # (R_k, S_k) bitmasks
cost = sum(rho(R,A) for R,A in terms); print('cost', cost)
n=120; m=8
Rs=[[k for k in range(m) if terms[k][0]>>t&1] for t in range(n)]
Ss=[[k for k in range(m) if terms[k][1]>>t&1] for t in range(n)]
OM=[ORTH[t] | (1<<t) for t in range(n)]
def rank_of(A,Bh):
    rows=[]
    for t in range(n):
        r=0
        for s in bits(OM[t]):
            if bin(A[t]&Bh[s]).count('1')&1: r|=1<<s
        rows.append(r)
    return rank2(rows)
rnd=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 1)
best=0
for rep in range(4):
    A=[sum(1<<k for k in Rs[t] if rnd.random()<.5) for t in range(n)]
    Bh=[sum(1<<k for k in Ss[t] if rnd.random()<.5) for t in range(n)]
    cur=rank_of(A,Bh)
    for it in range(6000):
        t=rnd.randrange(n)
        if rnd.random()<.5:
            k=rnd.choice(Rs[t]); A[t]^=1<<k; v=rank_of(A,Bh)
            if v>=cur: cur=v
            else: A[t]^=1<<k
        else:
            k=rnd.choice(Ss[t]); Bh[t]^=1<<k; v=rank_of(A,Bh)
            if v>=cur: cur=v
            else: Bh[t]^=1<<k
    print('rep',rep,'cert',cur, flush=True); best=max(best,cur)
print('best',best)
print('PASS_EXPNZ_CERTHC' if cost == 72 and best == 72 else 'FAIL_EXPNZ_CERTHC')   # the F2 packing certificate is tight (72 = cost) on the P120 frame decoder

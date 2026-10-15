"""Certificate: dim-2 register count R_2 >= 336 in the expSF OR-reach relaxation (=> >=112 helpers alive at dim 2).
Per port s, the 52 pairs (s,t) leave <p_s> only through dim-2 out-neighbours N of <p_s>; reach(N) <= ports(N^perp).
Charge m(N) (registers at N) to the ports whose line-edge into N is used (split m/2 for 2-port nodes).
Facts checked here: (a) 2-port nodes reach <=26; (b) the source's own L2 (u in K_s) reaches <=24 if 2-port, <=28 if 1-port;
(c) every dim-6 subspace holds <=20 ports, so a node with one out-edge (m=1) reaches <=20; (d) min charge per port >=3."""
import json, numpy as np, itertools
from fractions import Fraction as F
P=json.load(open('instA.json'))['ports']; cen=json.load(open('cen.json'))['cen']; n=112
nor=[1,7,25,31,121,127,97,103]; dot=lambda a,b:(a&b).bit_count()&1; isp={p:i for i,p in enumerate(P)}
T=lambda Q: sum(1 for p in P if all(dot(p,q)==0 for q in Q))
ok=lambda c,m: print(('PASS ' if c else 'FAIL ')+m) or c
a=max(T([P[s],P[q]]) for s in range(n) for q in range(n) if q!=s); ok(a==26,f'(a) max reach of a 2-port dim-2 node = {a}')
r2=0; r1=0; samecls=True
for s in range(n):
    for u in range(1,512):
        if u==P[s] or not all(dot(u,nor[q])==0 for q in cen[s]): continue
        two = (u in isp) or ((u^P[s]) in isp)
        if two:
            q=isp[u] if u in isp else isp[u^P[s]]; samecls&= cen[q]==cen[s]; r2=max(r2,T([P[s],u]))
        else: r1=max(r1,T([P[s],u]))
ok(samecls and r2==24 and r1==28,f'(b) source L2 in K_s: 2-port => same class, reach<={r2}; 1-port reach<={r1}')
M=np.zeros((512,2),np.uint64)
for q in range(512):
    m=sum(1<<i for i,p in enumerate(P) if dot(p,q)==0); M[q,0]=m&((1<<64)-1); M[q,1]=m>>64
mp6=0
for x in range(1,512):
    for y in range(x+1,512):
        z=np.arange(y+1,512); z=z[z!=(x^y)]
        if len(z): mp6=max(mp6,int(np.bitwise_count(M[x]&M[y]&M[z]).sum(1).max()))
ok(mp6==20,f'(c) max ports in any 6-dim subspace = {mp6}')
# (d) node types (reach cap, charge to s). generic out-neighbour / source's own node
gen=[(26,F(1)),(26,F(3,2)),(mp6,F(1)),(32,F(2)),(26,F(2))]   # 2-port shared m=2/3; 1-port or single-use m=1; m>=2
src=[(24,F(1)),(24,F(3,2)),(min(r1,mp6),F(1)),(r1,F(2)),(24,F(2))]
best=None
for sn in src:
    for k in range(0,4):
        for combo in itertools.combinations_with_replacement(gen,k):
            if sn[0]+sum(c[0] for c in combo)>=52:
                c=sn[1]+sum(c[1] for c in combo); best=c if best is None or c<best else best
ok(best==3,f'(d) minimum charge per port to reach all 52 tops = {best}  => R_2 >= {best*n}, helpers alive at dim 2 >= {best*n-224}')

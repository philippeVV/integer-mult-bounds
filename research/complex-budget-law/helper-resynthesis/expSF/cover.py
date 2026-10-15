"""Stage-1 necessary condition (hub-relaxed, zero non-constant helpers): x_s leaves the island <p_s> only along
chains that contain <p_s>: source s (next level G_s=<p_s,u_s>), the entry target e_s (next level <p_s,w_s>, w_s ⊥ p_e),
<=1 centre first-read per centre. Every reader t needs one of these next levels inside p_t^perp."""
import json
from collections import Counter
from itertools import combinations
I=json.load(open('instA.json')); P=I['ports']; al=json.load(open('alpha.json')); cen=json.load(open('cen.json'))['cen']
nor=[1,7,25,31,121,127,97,103]
dot=lambda a,b:(a&b).bit_count()&1
R={s:[int(t) for t in al if str(s) in al[t]] for s in range(112)}
orth={s:[t for t in range(112) if dot(P[s],P[t])==0] for s in range(112)}
print('readers == all orthogonal ports?',Counter(len(orth[s])-len(R[s]) for s in range(112)), 'self-orth ports',sum(dot(P[s],P[s])==0 for s in range(112)))
V=range(1,512)
def Kset(s): return [x for x in V if all(dot(x,nor[q])==0 for q in cen[s])]
def mincover(s, first):
    rp=[P[t] for t in R[s]]
    best=None
    for k in range(1,6):
        # greedy exact search: choose k vectors (first from 'first'), cover rp by perps
        def rec(rem,k,pool):
            if not rem: return []
            if k==0: return None
            for x in pool:
                if x in (P[s],): continue
                r2=[y for y in rem if dot(x,y)]
                if len(r2)==len(rem): continue
                sub=rec(r2,k-1,V)
                if sub is not None: return [x]+sub
            return None
        r=rec(rp,k,first)
        if r is not None: return k,r
for lab,first in (('u in K_s',None),('u free',list(V))):
    c=Counter(); ex={}
    for s in range(112):
        k,r=mincover(s, first if first else Kset(s)); c[k]+=1; ex[s]=r
    print(lab,'min #next-level vectors to cover readers(s):',sorted(c.items()))
# pairs: for each s, the set of w completing a cover with some u in K_s
W2={}
for s in range(112):
    rp=[P[t] for t in R[s]]; ws=set()
    for u in Kset(s):
        if u==P[s]: continue
        rem=[y for y in rp if dot(u,y)]
        for w in V:
            if w!=P[s] and all(dot(w,y)==0 for y in rem): ws.add(w)
    W2[s]=ws
print('#s with a 2-cover (u in K_s):',sum(1 for s in W2 if W2[s]))

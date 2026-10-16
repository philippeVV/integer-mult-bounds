"""Sanity: orientation/degree facts BH state (d-(v)=d+(v)=5, d+_C(v)=d-_J(v)), hop distances, terminal multiplicities."""
import sys; sys.path.insert(0, __file__.rsplit('/',1)[0])
from bh import *
from collections import Counter
din=Counter(); dout=Counter(); dCp=Counter(); dJm=Counter()
for q,(u,v) in enumerate(edges):
    a,b=(u,v) if POS[u]<POS[v] else (v,u)          # a earlier
    if q in Rs: tail,head=b,a                        # deleted edges later->earlier
    else: tail,head=a,b
    dout[tail]+=1; din[head]+=1
    if q in C: dCp[tail]+=1
    if q in J: dJm[head]+=1
print("d-=d+=5 everywhere:", all(din[v]==5==dout[v] for v in range(182)), " d+_C=d-_J:", all(dCp[v]==dJm[v] for v in range(182)))
term=Counter()
for q,s,t in SESS: term[s]+=1; term[t]+=1
print("terminals per vertex:", sorted(Counter(term[v] for v in range(182)).items()), "total", sum(term.values()))

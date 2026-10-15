import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""Total B+restore correction in original source variables, separated into post-cut target content and
orthogonal-source part: x_t(final) = sum_t' beta y_t' + sum_s alpha_ts x_s (helpers expanded to their A-content)."""
import sys, json; sys.path.insert(0,'../expBV')
from wordmodel import Word, basis, regs
from fractions import Fraction as F
from collections import Counter
W=Word('../../helper-bounds/sel1094.json.gz'); w=W.w; v=W.v; n=W.n; P=w['ports']
dot=lambda a,b:(a&b).bit_count()&1
def ap(e,rows):
    if e[0]=='add':
        t,s=e[1],e[2]; q=F(e[3],e[4])
        for k,x in list(rows[s].items()):
            a=rows[t].get(k,0)+q*x
            if a: rows[t][k]=a
            else: rows[t].pop(k,None)
    else:
        q=F(e[2],e[3]); rows[e[1]]={k:q*x for k,x in rows[e[1]].items()}
# A-content of every register (pre+comp+A), as original vars
A=[{k:F(1)} for k in range(n)]
for sec in('precondition','compensation','A'):
    for e in w[sec]: ap(e,A)
B=[{k:F(1)} for k in range(n)]
for e in w['B']: ap(e,B)
alpha={}; beta={}
for t in range(v):
    a=Counter(); b={}
    for k,x in B[v+t].items():
        if v<=k<2*v: b[k-v]=x; continue
        for j,y in A[k].items():
            if j<v: a[j]+=x*y
    alpha[t]={j:x for j,x in a.items() if x}; beta[t]=b
sz=Counter(len(a) for a in alpha.values()); print('alpha support sizes',sorted(sz.items()))
print('alpha orth', Counter(dot(P[s],P[t]) for t in alpha for s in alpha[t]))
print('beta sizes',Counter(len(b) for b in beta.values()))
rd=Counter(); 
for t in alpha:
    for s in alpha[t]: rd[s]+=1
print('reads per source',sorted(Counter(rd.values()).items()))
Kd=Counter()
for s in range(112):
    ts=[t for t in alpha if s in alpha[t]]
    K=basis([P[t] for t in ts]); Kd[9-len(K)]+=1
print("K'_s dims (meet of p_t^perp over readers)",sorted(Kd.items()))
json.dump({str(t):{str(s):str(x) for s,x in alpha[t].items()} for t in alpha},open('alpha.json','w'))
# span dims of alpha supports per target
print('span dim of alpha_t ports',sorted(Counter(len(basis([P[s] for s in alpha[t]])) for t in alpha).items()))

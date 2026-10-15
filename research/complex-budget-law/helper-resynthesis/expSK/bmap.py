import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
import sys, json; sys.path.insert(0,'../expBV')
from wordmodel import Word, basis, regs
from fractions import Fraction as F
from collections import Counter
W=Word('../../helper-bounds/sel1094.json.gz'); w=W.w; v=W.v; n=W.n; P=w['ports']
dot=lambda a,b:(a&b).bit_count()&1
# B map: start from identity at post-cut; apply B
rows=[{k:F(1)} for k in range(n)]
def ap(e,rows):
    if e[0]=='add':
        t,s=e[1],e[2]; q=F(e[3],e[4])
        for k,x in list(rows[s].items()):
            a=rows[t].get(k,0)+q*x
            if a: rows[t][k]=a
            else: rows[t].pop(k,None)
    else:
        q=F(e[2],e[3]); rows[e[1]]={k:q*x for k,x in rows[e[1]].items()}
for e in w['B']: ap(e,rows)
# frames at end of A (before B): replay frames up to start of B
tB=min(t for t,(s,e) in enumerate(W.events) if s=='B')
cur=list(W.start)
for t,(s,e) in enumerate(W.events[:tB]):
    if s in('precondition','compensation','cut') or e[0]=='warp': continue
    for q in regs(e): cur[q]=basis(e[-1])
# targets at cut level ()
tS=Counter(); tT=Counter(); tH=Counter(); orth=Counter(); hdim=Counter()
for t in range(v):
    r=rows[v+t]
    S_=[k for k in r if k<v]; T_=[k for k in r if v<=k<2*v and k!=v+t]; H_=[k for k in r if k>=2*v]
    tS[len(S_)]+=1; tT[len(T_)]+=1; tH[len(H_)]+=1
    for k in S_: orth[dot(P[k],P[t])]+=1
    for k in H_: hdim[len(cur[k])]+=1
print('B-map deps per target: sources',sorted(tS.items()),'\n targets',sorted(tT.items()),'\n helpers',sorted(tH.items()))
print('source deps p_s.p_t parity',orth,' helper post-A frame dims',sorted(hdim.items()))

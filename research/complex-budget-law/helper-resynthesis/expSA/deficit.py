import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""Island/direct-capacity accounting on the selected #1094 interface.
(1) First-touch lemma inputs: for every source s, meet of p_t^perp over its B-readers (alpha support) is <p_s>.
(2) If s climbs pre-cut to a frame containing <p_s,u> with u in K_s (A-side cap), only B-readers t with
    p_t . u = 0 can be reached by paths starting at s itself after its first non-island gate.
    deficit(s) = 52 - max_u #{t in readers(s): p_t.u = 0}; these pairs need an island entry (one per target)
    plus target-to-target relays if no helper is allowed to climb.
(3) current word: first gate on every source is at a dim-1 frame with a helper partner."""
import json, sys
from collections import Counter
sys.path.insert(0,'../expBV'); from wordmodel import Word, basis, regs
I=json.load(open('instA.json')); P=I['ports']
al=json.load(open('alpha.json')); cen=json.load(open('cen.json'))['cen']
nor=[1,7,25,31,121,127,97,103]
dot=lambda a,b:(a&b).bit_count()&1
readers={s:[int(t) for t in al if str(s) in al[t]] for s in range(112)}
meet1=all(len(basis([P[t] for t in readers[s]]))==8 and all(dot(P[s],P[t])==0 for t in readers[s]) for s in range(112))
print('(1) every source: 52 readers, all orthogonal, reader ports span p_s^perp (meet=<p_s>):',meet1, Counter(len(r) for r in readers.values()))
tot=0; best=Counter(); 
for s in range(112):
    K=[x for x in range(1,512) if all(dot(x,nor[q])==0 for q in cen[s])]
    m=max(sum(1 for t in readers[s] if dot(P[t],u)==0) for u in K if u!=P[s])
    best[m]+=1; tot+=52-m
print('(2) max direct readers after a 2-dim pre-cut climb inside K_s:',dict(best),' total deficit pairs',tot,' island entries available <=',112)
# with u ranging over ALL of F_2^9 (no A-side cap) for comparison
tot2=0
for s in range(112):
    m=max(sum(1 for t in readers[s] if dot(P[t],u)==0) for u in range(1,512) if u!=P[s]); tot2+=52-m
print('    (same, u unrestricted):',tot2)
W=Word('../../helper-bounds/sel1094.json.gz'); v=W.v
first={}
for t,(sec,e) in enumerate(W.events):
    if sec in('precondition','compensation','cut') or e[0]!='add': continue
    for q in regs(e):
        if q<v and q not in first:
            other=e[1] if e[2]==q else e[2]
            first[q]=(len(basis(e[-1])), 'H' if other>=2*v else ('S' if other<v else 'T'))
print('(3) current word first framed gate per source (dim, partner):',Counter(first.values()))

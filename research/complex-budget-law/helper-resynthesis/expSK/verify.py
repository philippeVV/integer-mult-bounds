import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""Independent OR-reach verification of a design file (relaxed hub model of expSF/expSH):
frames monotone (full flags, +1 dim per step), source starts at its port line with dim-2 frame inside K_s,
target ends at p_t^perp, content born only in source registers, merged at every frame value (cluster rule),
every target top holds the 52 alpha_t sources, every centre frame U_q holds its 56 sources."""
import json, sys, math
from collections import Counter, defaultdict
sys.path.insert(0,'../expBV'); from wordmodel import basis
D=json.load(open(sys.argv[1])); P=json.load(open('instA.json'))['ports']; I=json.load(open('instA.json'))
cen=json.load(open('cen.json'))['cen']; al=json.load(open('alpha.json')); n=112
nor=[1,7,25,31,121,127,97,103]; dot=lambda a,b:(a&b).bit_count()&1
perp=lambda B: basis([x for x in range(1,512) if all(dot(x,b)==0 for b in B)])
contained=lambda a,b: basis(list(a)+list(b))==tuple(b)
ok=lambda c,m: print(('PASS ' if c else 'FAIL ')+m) or c
R2='--r2' in sys.argv  # audit rule: non-targets never read targets
print('R2 (targets only absorb):',R2)
regs=[[tuple(N) for N in p] for p in D['pieces']]; role=[tuple(r) if r!='H' else 'H' for r in D['role']]
ok(all(contained(a,b) and len(b)==len(a)+1 for ch in regs for a,b in zip(ch,ch[1:])),'all chains are monotone full flags')
K=[perp([nor[q] for q in cen[s]]) for s in range(n)]
src={r[1]:i for i,r in enumerate(role) if r!='H' and r[0]=='S'}; tgt={r[1]:i for i,r in enumerate(role) if r!='H' and r[0]=='T'}
ok(sorted(src)==list(range(n)) and sorted(tgt)==list(range(n)),'112 sources and 112 targets, one register each')
ok(all(regs[i][0]==(P[s],) for s,i in src.items()),'each source register starts at its port line <p_s>')
ok(all(len(regs[i])>1 and contained(regs[i][1],K[s]) for s,i in src.items()),'each source dim-2 frame lies in K_s (centre-hyperplane meet)')
ok(all(regs[i][-1]==perp([P[t]]) for t,i in tgt.items()),'each target register ends at p_t^perp')
ok(all(contained(N,perp([P[t]])) for t,i in tgt.items() for N in regs[i]),'target chains stay inside p_t^perp')
# content propagation by dimension
cont=[set([r[1]]) if (r!='H' and r[0]=='S') else set() for r in role]
at=defaultdict(list)
for i,ch in enumerate(regs):
    for k,N in enumerate(ch): at[(len(N),N)].append(i)
nodecont={}
for (dm,N) in sorted(at,key=lambda x:x[0]):
    c=set(); ct=set()
    for i in at[(dm,N)]:
        if R2 and role[i]!='H' and role[i][0]=='T': ct|=cont[i]
        else: c|=cont[i]
    nodecont[N]=c|ct
    for i in at[(dm,N)]: cont[i]=set(c|ct) if (role[i]!='H' and role[i][0]=='T') or not R2 else set(c)
ok(all(all(dot(P[s],x)==0 for x in perp(N) ) for N,c in nodecont.items() for s in c),'framed content: every carried port lies in the frame')
good=0
for t in range(n):
    need={int(s) for s in al[str(t)]}; got=nodecont.get(perp([P[t]]),set())
    good+=len(need&got)
ok(good==5824,f'target tops deliver {good}/5824 (s,t) alpha pairs')
U=[basis(c['U']) for c in I['centres']]
cg=0
for q in range(8):
    need={s for s in range(n) if dot(P[s],nor[q])==0}; got=nodecont.get(U[q],set()); cg+=len(need&got)
    assert perp([nor[q]])==U[q]
ok(cg==448,f'centre frames U_q receive {cg}/448 (s,q) pairs via non-constant registers')
f=lambda w: w*math.log(73/w)
hw=Counter(len(regs[i][-1])-len(regs[i][0]) for i,r in enumerate(role) if r=='H')
print('registers',len(regs),'helpers',sum(hw.values()),'inventory',dict(sorted(hw.items())),'log-mass %.1f'%sum(k*f(w) for w,k in hw.items()))
print('registers per dim',dict(sorted(Counter(len(N) for ch in regs for N in ch).items())))

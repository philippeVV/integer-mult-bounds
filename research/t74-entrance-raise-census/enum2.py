import os
import sys,pickle; sys.path.insert(0,'.')
from classes import *
WD=os.environ.get('PKG','gmx-t74-borrow43-handoff')+'/new/borrow/words/'
nm=sys.argv[1]
W=Word(load(WD+'BORROW-WORD-43%s.json.gz'%nm)); C=Classes(W)
print('classes',len(C.members))
seedcls={}
for h in C.ent:
    s=C.seq.get(h)
    if s and s[0][0]=='C': seedcls.setdefault(s[0][1],[]).append(h)
print('seed classes',len(seedcls))
out=[]
for c0,hs in seedcls.items():
    F=C.frame[c0]; S=span_set(F); seen=set()
    for x in range(1,512):
        if x in S: continue
        key=basis(list(F)+[x])
        if key in seen: continue
        seen.add(key)
        L,flags=C.closure([c0],x)
        if flags: out.append(dict(c0=c0,x=x,flags=sorted(flags),nL=len(L))); continue
        ev=C.evaluate(L,x)
        if ev is None: out.append(dict(c0=c0,x=x,flags=['warpseq'],nL=len(L))); continue
        dp,raised=ev
        out.append(dict(c0=c0,x=x,flags=[],nL=len(L),dp=dp,raised=raised,L=sorted(L)))
pickle.dump(dict(out=out,seedcls=seedcls),open('enum2-%s.pkl'%nm,'wb'))
from collections import Counter
print(Counter(tuple(o['flags']) for o in out))
ok=[o for o in out if not o['flags']]
for o in sorted(ok,key=lambda o:-len(o['raised']))[:40]:
    print(o['c0'],C.frame[o['c0']],o['x'],'L',o['nL'],'raised',len(o['raised']),'dp',sorted(o['dp'].items()))

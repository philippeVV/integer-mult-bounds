import os
import sys,pickle; sys.path.insert(0,'.')
from classes import *
WD=os.environ.get('PKG','gmx-t74-borrow43-handoff')+'/new/borrow/words/'
W=Word(load(WD+'BORROW-WORD-43extra.json.gz')); C=Classes(W)
d=pickle.load(open('enum2-extra.pkl','rb'))
from collections import Counter
import statistics
print('nL dist', Counter(min(o['nL']//5*5,100) for o in d['out']))
# per seed class: minimum nL, and end-cause
for c0,hs in list(d['seedcls'].items())[:52]:
    os_=[o for o in d['out'] if o['c0']==c0]
    best=min(os_,key=lambda o:o['nL'])
    x=best['x']
    L,fl=C.closure([c0],x)
    # causes
    causes=[]
    for c in L:
        for r in C.regs[c]:
            s=C.seq[r]; i=[j for j,it in enumerate(s) if it==('C',c)][0]
            if i+1==len(s) and not inspan(x,W.end[r]): causes.append(('T' if 112<=r<224 else 'H',r,len(W.end[r])))
            if i+1<len(s) and s[i+1][0]=='W': causes.append(('W',r))
    print(c0,C.frame[c0],'helpers',sorted(hs),'minL',best['nL'],'x',x,best['flags'],causes[:6])

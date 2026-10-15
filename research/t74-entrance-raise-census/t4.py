import os
import sys,pickle; sys.path.insert(0,'.')
from frames import *
st=pickle.load(open('states.pkl','rb'))['extra']
A=Word(st['handoff']); B=Word(st['entrance'])
raised=[int(k) for k in B.w['start_extra'] if B.w['start_extra'][k]!=A.w['start_extra'].get(k)]
print('raised',raised)
newst={h:B.start[h] for h in raised}
for h in raised:
    fr,fl,t=A.lift({h:B.start[h]}); print(h,A.start[h],B.start[h],sorted(fl)[:5],len(fr))
fr,fl,t=A.lift(newst); print('joint',sorted(fl)[:8],len(fr))
# compare entrance frames vs handoff frames: which events changed
ch=[k for k in range(len(A.ev)) if A.ev[k][2]!=B.ev[k][2]]
print('changed events',len(ch))
sup=sum(1 for k in ch if A.ev[k][0]=='gate' and contained(A.ev[k][2],B.ev[k][2]))
print('changed that are supersets',sup)
for k in ch[:30]: print(k,A.ev[k][1],A.ev[k][2],'->',B.ev[k][2])

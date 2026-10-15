import os
import sys,pickle; sys.path.insert(0,'.')
from frames import *
WD=os.environ.get('PKG','gmx-t74-borrow43-handoff')+'/new/borrow/words/'
res={}
for nm in ('base','extra'):
    W=Word(load(WD+'BORROW-WORD-43%s.json.gz'%nm))
    out=[]
    for h in sorted(map(int,W.w['start_extra'])):
        s=W.start[h]; seen=set()
        for x in range(1,512):
            S2=basis(list(s)+[x])
            if len(S2)==len(s) or S2 in seen: continue
            seen.add(S2)
            frames,flags,touched=W.lift({h:S2})
            dp=Counter()
            for r in touched:
                for d in W.climbs(r): dp[d]-=1
                for d in W.climbs(r,frames,S2 if r==h else None): dp[d]+=1
            dp=+dp+Counter(); dpn={k:v for k,v in (Counter(dp)).items()}
            # full signed delta
            sd=Counter()
            for r in touched:
                for d in W.climbs(r): sd[d]-=1
                for d in W.climbs(r,frames,S2 if r==h else None): sd[d]+=1
            sd={k:v for k,v in sd.items() if v}
            out.append(dict(h=h,S=S2,flags=sorted(flags),dp=sd,nch=len(frames),touched=len(touched),w=len(W.end[h])-len(s)))
    res[nm]=out
    print(nm,len(out))
pickle.dump(res,open('enum1.pkl','wb'))

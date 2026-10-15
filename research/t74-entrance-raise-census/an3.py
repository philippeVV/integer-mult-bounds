import os
import pickle,sys; sys.path.insert(0,'.')
from frames import *
from collections import Counter,defaultdict
res=pickle.load(open('enum1.pkl','rb'))
W=Word(load(os.environ.get('PKG','gmx-t74-borrow43-handoff')+'/new/borrow/words/BORROW-WORD-43extra.json.gz'))
by=defaultdict(list)
for c in res['extra']: by[c['h']].append(c)
shape=lambda h:(len(W.start[h]),len(W.end[h])-len(W.start[h]),(W.climbs(h) or [None])[0])
tab=Counter(); feas=[]
for h,cs in sorted(by.items()):
    kinds=Counter()
    for c in cs:
        f=set(x.split(':')[0] for x in c['flags'])
        if not f: kinds['ok']+=1; feas.append(c)
        elif 'warp' in f: kinds['warp+end']+=1
        else:
            # end causes: targets vs cold helpers
            ends=[int(x.split(':')[1]) for x in c['flags'] if x.startswith('end:')]
            t=any(112<=r<224 for r in ends); hh=any(r>=224 for r in ends); xs=any(r<112 for r in ends)
            kinds['end-'+('T' if t else '')+('H' if hh else '')+('X' if xs else '')]+=1
    first=W.byreg[h][0]; partners=[r for r in W.ev[first][1] if r!=h]
    ptype=['T' if 112<=r<224 else ('X' if r<112 else ('E' if str(r) in W.w['start_extra'] else 'H')) for r in partners]
    tab[(shape(h),ptype[0] if ptype else None,tuple(sorted(kinds)))]+=1
for k,v in sorted(tab.items(),key=str): print(v,k)
print('feasible single raises')
for c in feas: print(c['h'],c['S'],c['dp'],'events',c['nch'])
pickle.dump(feas,open('feas.pkl','wb'))

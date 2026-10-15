import os
import pickle,sys; sys.path.insert(0,'.')
from frames import *
from lp import *
W=Word(load(os.environ.get('PKG','gmx-t74-borrow43-handoff')+'/new/borrow/words/BORROW-WORD-43extra.json.gz'))
feas=pickle.load(open('feas.pkl','rb'))
r0=root(INV0,P0,lo=0.00165,hi=0.00175,it=50)
def price(dp,wch):
    inv=dict(INV0)
    for a,b in wch: inv[a]-=1; inv[b]=inv.get(b,0)+1
    p=Counter(P0); p.update(dp)
    return root(inv,dict(p),lo=0.00165,hi=0.00175,it=50)-r0
for c in feas:
    w0=len(W.end[c['h']])-len(W.start[c['h']])
    if w0-1 not in LP.WID: print(c['h'],'width',w0,'->',w0-1,'not in LP widths'); continue
    print(c['h'],c['dp'],'w',w0,'->',w0-1,'dkappa(LP) %.3e'%price(c['dp'],[(w0,w0-1)]))

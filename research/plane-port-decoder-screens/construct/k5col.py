import json, sys
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType
K5=json.load(open('k5.json')); n=112; m=int(sys.argv[1]); p=int(sys.argv[2]) if len(sys.argv)>2 else 2
# vars: x[t][j] colour j (one-hot), d[t] drop
x=[[t*p+j+1 for j in range(p)] for t in range(n)]; d=[n*p+t+1 for t in range(n)]
s=Cadical153(); top=n*p+n
for t in range(n):
    s.add_clause(x[t]+[d[t]])
for k in K5:
    for j in range(p):
        s.add_clause([-x[t][j] for t in k])
s.add_clause([x[0][0]])
cnf=CardEnc.atmost(lits=d,bound=m,top_id=top,encoding=EncType.seqcounter)
for cl in cnf.clauses: s.add_clause(cl)
print('p',p,'m',m,s.solve())

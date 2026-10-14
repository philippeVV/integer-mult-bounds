"""SAT: 2-class tensor on Q = ports minus dropped set P, |P|<=m. args: m [restrict8]"""
import sys, json; sys.path.insert(0,'.')
from csp import *
import pycryptosat
from pysat.card import CardEnc, EncType
m=int(sys.argv[1]); restrict=len(sys.argv)>2 and sys.argv[2]=='r8'
nv=[0]
def new():
    nv[0]+=1; return nv[0]
s=pycryptosat.Solver(threads=4)
c=[new() for _ in range(n)]
r=[[new() for _ in range(9)] for _ in range(n)]
dr=[new() for _ in range(n)]
Gd=[(G[i]>>i)&1 for i in range(9)]
def xor(vs,rhs): s.add_xor_clause(vs,rhs)
Cf={i:[ip(L[ports[i]],L[f]) for f in fr] for i in range(n)}
for t in range(n):
    Gu=[gform(1<<i,U[t],G) for i in range(9)]
    xor([r[t][i] for i in range(9) if Gd[i]^Gu[i]],True)
    if restrict:
        for k in range(8):
            if Cf[t][k]:
                Gf=[gform(1<<i,labs[fr[k]],G) for i in range(9)]
                xor([r[t][i] for i in range(9) if Gf[i]],False)
g=[[new() for _ in range(9)] for _ in range(n)]
for t in range(n):
    for i in range(9):
        xor([g[t][i]]+[r[t][j] for j in range(9) if G[i]>>j&1],False)
s.add_clause([-c[0]])
for a in range(n):
    for b in odd[a]:
        if b<a: continue
        e=new()
        xor([e,c[a],c[b]],True)
        act=[-e,dr[a],dr[b]]
        for (x,y) in ((a,b),(b,a)):
            yv=new(); Gu=[gform(1<<i,U[x],G) for i in range(9)]
            xor([yv]+[r[y][i] for i in range(9) if Gu[i]],False)
            s.add_clause(act+[-yv])
        ps=[]
        for i in range(9):
            p=new(); ps.append(p)
            s.add_clause([-p,r[a][i]]); s.add_clause([-p,g[b][i]]); s.add_clause([p,-r[a][i],-g[b][i]])
        y3=new(); xor([y3]+ps,False); s.add_clause(act+[-y3])
cnf=CardEnc.atmost(lits=dr,bound=m,top_id=nv[0],encoding=EncType.seqcounter)
for cl in cnf.clauses: s.add_clause(cl)
sat,sol=s.solve()
print('m',m,'restrict',restrict,'SAT' if sat else 'UNSAT')
if sat:
    out={'cls':[int(sol[c[t]]) for t in range(n)],'R':[sum(int(sol[r[t][i]])<<i for i in range(9)) for t in range(n)],
         'drop':[int(sol[dr[t]]) for t in range(n)]}
    json.dump(out,open(f'sat3_m{m}_{int(restrict)}.json','w'))
    print('dropped',sum(out['drop']))

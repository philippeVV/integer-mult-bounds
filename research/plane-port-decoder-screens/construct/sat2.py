"""SAT: 2-class partition + plane r_t with all same-class odd pairs plane-orthogonal."""
import sys, json; sys.path.insert(0,'.')
from csp import *
import pycryptosat
nv=[0]
def new():
    nv[0]+=1; return nv[0]
s=pycryptosat.Solver(threads=4)
c=[new() for _ in range(n)]
r=[[new() for _ in range(9)] for _ in range(n)]
bits=lambda v:[i for i in range(9) if v>>i&1]
Gd=[(G[i]>>i)&1 for i in range(9)]
def xor(vs,rhs): s.add_xor_clause(vs,rhs)
for t in range(n):
    # r.r + u.r = 1 : sum_i (Gd_i + (G u)_i) r_i
    Gu=[gform(1<<i,U[t],G) for i in range(9)]
    xor([r[t][i] for i in range(9) if Gd[i]^Gu[i]],True)
g=[[new() for _ in range(9)] for _ in range(n)]
for t in range(n):
    for i in range(9):  # g_t[i] = sum_j G_ij r_t[j]
        xor([g[t][i]]+[r[t][j] for j in range(9) if G[i]>>j&1],False)
s.add_clause([-c[0]])
for a in range(n):
    for b in odd[a]:
        if b<a: continue
        e=new()  # e -> c_a==c_b ;  c_a==c_b -> e
        xor([e,c[a],c[b]],True)  # e = 1 xor c_a xor c_b
        for (x,y) in ((a,b),(b,a)):
            yv=new(); Gu=[gform(1<<i,U[x],G) for i in range(9)]
            xor([yv]+[r[y][i] for i in range(9) if Gu[i]],False)
            s.add_clause([-e,-yv])
        ps=[]
        for i in range(9):
            p=new(); ps.append(p)
            s.add_clause([-p,r[a][i]]); s.add_clause([-p,g[b][i]]); s.add_clause([p,-r[a][i],-g[b][i]])
        y3=new(); xor([y3]+ps,False); s.add_clause([-e,-y3])
extra=sys.argv[1] if len(sys.argv)>1 else ''
if extra=='balanced':
    pass
sat,sol=s.solve()
print('SAT' if sat else 'UNSAT', 'vars',nv[0])
if sat:
    cls=[int(sol[c[t]]) for t in range(n)]
    R=[sum(int(sol[r[t][i]])<<i for i in range(9)) for t in range(n)]
    json.dump({'cls':cls,'R':R},open('sat2_sol.json','w'))
    print(sum(cls))

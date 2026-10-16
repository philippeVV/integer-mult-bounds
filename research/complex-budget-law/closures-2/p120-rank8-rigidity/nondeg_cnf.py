"""Independent cross-check of the nondegeneracy lemma (README Step 1): Z1-only Boolean system, plain CNF, solved by
CaDiCaL/Glucose (PySAT).  By default it rebuilds the CNF, checks it is identical to the stored nondeg.cnf.gz and that both
solvers return UNSAT.  With --write it overwrites nondeg.cnf.gz and the Glucose DRUP proof nondeg.drup.gz.
Run: uv run -q --no-project --python 3.12 --with numpy --with python-sat python -I -B nondeg_cnf.py [--write]"""
import gzip, itertools, os, sys, tempfile, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__))
WRITE = '--write' in sys.argv
from pysat.formula import CNF, IDPool
from pysat.solvers import Solver
from e8 import R, G, n
ORTH=(G==0)
Q=[]
for a in range(n):
    if all(ORTH[a,b] for b in Q): Q.append(a)
    if len(Q)==8: break
inQ=set(Q); C=(R@R[Q].T)//4
S={a:[j for j in range(8) if C[a,j]!=0] for a in range(n)}
pool=IDPool(); cnf=CNF()
TRUE=pool.id('T'); cnf.append([TRUE])
def X(a,j): return TRUE if a in inQ else pool.id(('x',a,j))
def Y(a,j): return TRUE if a in inQ else pool.id(('y',a,j))
def AND(u,v):
    p=pool.id(('p',u,v)); cnf.extend([[-p,u],[-p,v],[p,-u,-v]]); return p
for a in range(n):
    for b in range(n):
        com=[s for s in S[a] if s in S[b]]
        if not com: continue
        if a==b: cnf.append([AND(X(a,s),Y(a,s)) for s in com]); continue
        if not ORTH[a,b]: continue
        ps=[AND(X(a,s),Y(b,s)) for s in com]
        for k,p in enumerate(ps): cnf.append([-p]+[q for l,q in enumerate(ps) if l!=k])   # count != 1
cnf.append([-X(a,j) for a in range(n) if a not in inQ for j in S[a]]+[-Y(a,j) for a in range(n) if a not in inQ for j in S[a]])
tmp = os.path.join(tempfile.mkdtemp(), 'nondeg.cnf'); cnf.to_file(tmp); text = open(tmp).read()
print('vars',pool.top,'clauses',len(cnf.clauses))
stored = os.path.join(HERE, 'nondeg.cnf.gz')
if WRITE:
    with gzip.open(stored, 'wt') as g: g.write(text)
same = gzip.open(stored, 'rt').read() == text
print('identical to stored nondeg.cnf.gz:', same)
unsat = 0
for name in ('cadical153','glucose4'):
    try:
        with Solver(name=name,bootstrap_with=cnf.clauses,with_proof=(name=='glucose4')) as s:
            r=s.solve(); print(name,'SAT' if r else 'UNSAT'); unsat += not r
            if not r and name=='glucose4' and WRITE:
                pr=s.get_proof()
                with gzip.open(os.path.join(HERE, 'nondeg.drup.gz'), 'wt') as g: g.write('\n'.join(pr)+'\n')
                print('proof lines',len(pr))
    except Exception as e: print(name,'error',e)
print('PASS_NONDEG_CNF' if same and unsat == 2 else 'FAIL_NONDEG_CNF')

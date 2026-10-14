"""E8 root lines, P112 ports, binary labels u_t in F2^9 (bitmask ints). Exact."""
import itertools, random
from fractions import Fraction as Q

def roots():
    R = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = [0]*8; v[i] = 2*si; v[j] = 2*sj; R.append(tuple(v))
    for s in itertools.product((1, -1), repeat=8):
        if s.count(-1) % 2 == 0: R.append(tuple(s))
    return R  # doubled coordinates: <a,b> = dot/4

def ip(a, b): return sum(x*y for x, y in zip(a, b)) // 4

def lines():
    L = []
    for r in roots():
        if tuple(-x for x in r) not in L: L.append(r)
    assert len(L) == 120
    return L

def dot(a, b): return bin(a & b).count('1') & 1

def rank2(vs):
    basis = []
    for v in vs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return len(basis)

def span_basis(vs):
    basis = []
    for v in vs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return basis

def form_rank(U, V):
    """rank of pairing matrix between bases of span U and span V."""
    bu, bv = span_basis(U), span_basis(V)
    rows = [sum(dot(a, b) << j for j, b in enumerate(bv)) for a in bu]
    return rank2(rows)


def solve_coords(B, v):
    """exact coords of v in basis B (list of 8 vectors) via Fractions."""
    n=8
    A=[[Q(B[j][i]) for j in range(n)]+[Q(v[i])] for i in range(n)]
    for c in range(n):
        p=next(r for r in range(c,n) if A[r][c]!=0); A[c],A[p]=A[p],A[c]
        for r in range(n):
            if r!=c and A[r][c]!=0:
                f=A[r][c]/A[c][c]; A[r]=[x-f*y for x,y in zip(A[r],A[c])]
    return [A[i][n]/A[i][i] for i in range(n)]

def simple_roots(L):
    # pick 8 roots forming a Z-basis: use standard E8 simple roots (doubled coords)
    s=[(1,-1,-1,-1,-1,-1,-1,1)]
    s+=[tuple(2 if k==0 else 2 if k==1 else 0 for k in range(8))]
    for i in range(6):
        v=[0]*8; v[i+1]=2; v[i]=-2; s.append(tuple(v))
    # sign: e_{i+1}-e_i ; and e0+e1
    return s

GM=None
def setup():
    """returns lines, labels (9-bit ints), Gram matrix GM (list of 9 row bitmasks), gform(a,b)."""
    L=lines(); B=simple_roots(L)
    # check unimodular Z-basis: all roots integer coords
    X=[]
    for r in L:
        c=solve_coords(B,r); assert all(x.denominator==1 for x in c), c
        X.append(sum((int(x)&1)<<i for i,x in enumerate(c)))
    G2=[[ip(B[i],B[j])&1 for j in range(8)] for i in range(8)]
    rows=[sum(G2[i][j]<<j for j in range(8)) for i in range(8)]
    rows.append(1<<8)
    labs=[x|(1<<8) for x in X]
    return L,labs,rows

def gform(a,b,rows):
    s=0
    i=0
    while a:
        if a&1: s^=bin(rows[i]&b).count('1')&1
        a>>=1; i+=1
    return s

def span_b(vs):
    basis=[]
    for v in vs:
        for b in basis: v=min(v,v^b)
        if v: basis.append(v)
    return basis

def frank(U,V,G):
    bu,bv=span_b(U),span_b(V)
    rows=[sum(gform(a,b,G)<<j for j,b in enumerate(bv)) for a in bu]
    return rank2(rows)

FRAME=[(2,2,0,0,0,0,0,0),(2,-2,0,0,0,0,0,0),(0,0,2,2,0,0,0,0),(0,0,2,-2,0,0,0,0),
       (0,0,0,0,2,2,0,0),(0,0,0,0,2,-2,0,0),(0,0,0,0,0,0,2,2),(0,0,0,0,0,0,2,-2)]

def p112():
    L,labs,G=setup()
    fr=[L.index(f) if f in L else L.index(tuple(-x for x in f)) for f in FRAME]
    ports=[i for i in range(120) if i not in fr]
    return L,labs,G,fr,ports

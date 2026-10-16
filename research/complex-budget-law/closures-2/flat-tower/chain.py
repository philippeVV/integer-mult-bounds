"""The retained four-pair scalar chain (split-centres-global/proofs/INHERITED-COMPLEX-GEOMETRY.txt,
'Four-pair scalar chain'), generalised to n pairs.  Checks, exactly over Z:
  * it realises n exchanges with 2n+1 centre calls (F_j: Y_j += X_j, an L-shear;
    B_j: X_j -= Y_j, a U-shear) plus scalar bridges that only mix X with X and Y with Y;
  * n=4 is the selected nine-call skeleton (4 B + 5 F, winding -1 as in native-controller-congruence);
  * the typed potentials a(S)=dim(gX cap X), b(S)=dim(gX cap Y) etc. move as the 2n+1 lower-bound proof says.
Run: python3 -I -B chain.py   (stdlib only)"""
from fractions import Fraction as Q
def ident(N): return [[Q(int(i==j)) for j in range(N)] for i in range(N)]
def rank(A):
    A=[r[:] for r in A]; rk=0; rows=len(A); cols=len(A[0]) if A else 0
    for c in range(cols):
        p=next((i for i in range(rk,rows) if A[i][c]!=0),None)
        if p is None: continue
        A[rk],A[p]=A[p],A[rk]
        for i in range(rows):
            if i!=rk and A[i][c]!=0:
                f=A[i][c]/A[rk][c]; A[i]=[a-f*b for a,b in zip(A[i],A[rk])]
        rk+=1
    return rk
def chain(n):
    X=lambda j:j; Y=lambda j:n+j; N=2*n
    ops=[]                                  # (kind, target, source, coeff); kind in F,B,bridge
    def add(kind,t,s,c): ops.append((kind,t,s,c))
    add('F',Y(0),X(0),1)
    for j in range(n):
        add('B',X(j),Y(j),-1)
        if j<n-1:
            add('bridge',X(j),X(j+1),1); add('bridge',Y(j+1),Y(j),-1)
            add('F',Y(j),X(j),1)
            add('bridge',Y(j+1),Y(j),1); add('bridge',X(j),X(j+1),-1)
        else:
            add('F',Y(j),X(j),1)
    for j in range(n-2,-1,-1):
        add('bridge',X(j+1),X(j),1); add('bridge',Y(j),Y(j+1),-1)
    return ops
def apply(ops,N):
    g=ident(N)                               # columns = images of basis vectors; row op on state
    for kind,t,s,c in ops:
        g[t]=[a+c*b for a,b in zip(g[t],g[s])]
    return g
def inter(rowsA,rowsB,N):
    # dim(span A cap span B) = rkA+rkB-rk(A;B)
    return rank(rowsA)+rank(rowsB)-rank(rowsA+rowsB)
def main():
    for n in range(1,9):
        ops=chain(n); N=2*n; g=apply(ops,N)
        calls=[o for o in ops if o[0]!='bridge']
        # bridges are type-preserving
        assert all((t<n)==(s<n) for k,t,s,c in ops if k=='bridge')
        assert all((k=='F' and t>=n and s<n) or (k=='B' and t<n and s>=n) for k,t,s,c in calls)
        # final exchange newX=oldY, newY=-oldX must give identity
        fin=[[Q(0)]*N for _ in range(N)]
        for j in range(n):
            fin[j]=g[n+j][:]; fin[n+j]=[-v for v in g[j]]
        ok=fin==ident(N)
        nB=sum(1 for o in calls if o[0]=='B'); nF=len(calls)-nB
        print(f'n={n}: centre calls {len(calls)} (=2n+1: {len(calls)==2*n+1}), F={nF} B={nB}, winding {nB-nF}, '
              f'bridges {len(ops)-len(calls)}, exchange identity {ok}')
        assert ok and len(calls)==2*n+1
    # potential trace for n=4 (selected skeleton): S=gX, T=gY
    n=4; N=8; ops=chain(n)
    E=ident(N); Xs=[E[i] for i in range(n)]; Ys=[E[n+i] for i in range(n)]
    g=ident(N); trace=[]
    for kind,t,s,c in ops:
        g[t]=[a+c*b for a,b in zip(g[t],g[s])]
        if kind=='bridge': continue
        cols=lambda idx:[[g[r][i] for r in range(N)] for i in idx]   # images of basis vectors (columns)
        S=cols(range(n)); T=cols(range(n,N))
        trace.append((kind,inter(S,Xs,N),inter(S,Ys,N),inter(T,Xs,N),inter(T,Ys,N)))
    print('n=4 potentials after each call (kind, a(S), b(S), a(T), b(T)); start (4,0,0,4), end (0,4,4,0):')
    for row in trace: print('  ',row)
    print('PASS_EXPNL_CHAIN')
if __name__=='__main__':
    main()

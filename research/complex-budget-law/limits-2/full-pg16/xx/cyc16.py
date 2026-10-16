"""Singer cyclic orders on full PG(2,16) with k_P=8 (points), k_L=9 (lines).
Point i ~ line (i+d) mod N, d in D.  Order: P_i at 2(m i mod N), L_l at 2((m l - c) mod N)+1 (cut before P_0).
Class orientation: displacement y_d=(m d - c) mod N; class out of P (P->L) is back iff wraps: count y_d; class into P: N-y_d.
F_class(m,c)= min over 9-out/8-in class choice = sum 9 smallest y + sum 8 largest (N-y)."""
import numpy as np, sys
from math import gcd
from singer import singer
N,D,_=singer(); D=np.array(D)
res=[]
cs=np.arange(N)
for m in range(1,N):
    if gcd(m,N)!=1: continue
    Y=np.sort((m*D[None,:]-cs[:,None])%N,axis=1)
    F=Y[:,:9].sum(1)+(N-Y[:,9:]).sum(1)
    c=int(F.argmin()); res.append((int(F[c]),m,c))
res.sort(); print('best class-F', res[:12]); print('distinct top values',sorted(set(r[0] for r in res))[:10])

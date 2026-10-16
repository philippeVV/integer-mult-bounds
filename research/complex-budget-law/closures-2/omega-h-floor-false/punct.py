"""Puncture P120 under a frame-gauge G/2 decoder; track cost, omega_req(Q), rank_F2 B[Q,Q]."""
import random, sys
from e8 import *
def omega(L):
    """clique number of the required graph (pairwise orthogonal lines, B=1) on line set L"""
    best=[0]
    def ext(sz,cand):
        if sz>best[0]: best[0]=sz
        while cand:
            if sz+bin(cand).count('1')<=best[0]: return
            c=cand.bit_length()-1; cand&=~(1<<c)
            ext(sz+1, cand & ORTH[c])
    ext(0,L); return best[0]
def rkQ(L): return rank2([BROW[a]&L for a in bits(L)])
def main():
    fr=frames(); M=[ROOT[a] for a in fr[0]]; T8=gauge_terms(M)
    def cost(L): return sum(rho(R&L,A&L) for R,A in T8)
    rnd=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 0)
    L=FULL; print('P120', cost(L), omega(L), rkQ(L))
    while bin(L).count('1')>40:
        best=None
        for a in bits(L):
            key=(cost(L&~(1<<a)), rnd.random())
            if best is None or key<best[0]: best=(key,a)
        L&=~(1<<best[1]); n=bin(L).count('1'); c=best[0][0]; w=omega(L); h=rkQ(L)
        print(n, c, w, h, 'min(V,wh)=',min(n,w*h), 'VIOL' if c<min(n,w*h) else '')
if __name__=="__main__": main()

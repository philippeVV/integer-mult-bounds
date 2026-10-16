import itertools, numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from geom import singer
N, D = singer(17); n=18; k=9; D=np.array(D); SUB=np.array(list(itertools.combinations(range(n),k)))
G = np.exp(2j*np.pi*np.outer(np.arange(k), D)/N)
c=np.ones(n)
for _ in range(2000):
    Q,_=np.linalg.qr((G*c).conj().T); l=(abs(Q)**2).sum(1); c*=np.sqrt(.5/l)
Q,_=np.linalg.qr((G*c).conj().T); print('freq{0..8} flattened minlogPl', np.linalg.slogdet(Q[SUB])[1].min(), 'median', np.median(np.linalg.slogdet(Q[SUB])[1]))
rng=np.random.default_rng(0); G=rng.normal(size=(k,n))+1j*rng.normal(size=(k,n)); Q,_=np.linalg.qr(G.conj().T); print('random', np.linalg.slogdet(Q[SUB])[1].min())

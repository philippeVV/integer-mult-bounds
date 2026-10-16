import numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from xg_core import *
from fractions import Fraction as F
from math import comb
# free (arcsine) values of the refined bound for J = 1..7  -> compare 1/g with g = 2J+2
for J in range(1, 8):
    m = np.array([1.0] + [comb(2*k, k) / 4**k / 2 for k in range(1, J + 1)])
    print('J', J, 'girth', 2*J+2, 'moment-only M =', round(markov(m, J), 6), ' refined =', round(markov_ref(m, J), 6), ' 1/g =', round(1/(2*J+2), 6))
V, nP, E, g = load('PG3'); Q = rand_codes(V, np.random.default_rng(5)); x = np.random.default_rng(1).standard_normal(E)*0.3
th, gr = gauge_value_grad(V, nP, E, Q, 2, x); h=1e-6
fd = np.array([(theta(V,nP,E,Q,2,x+h*np.eye(E)[i]) - theta(V,nP,E,Q,2,x-h*np.eye(E)[i]))/(2*h) for i in range(E)])
print('refined grad err', np.max(abs(fd-gr)), 'fourier PG3 refined', theta(V,nP,E,fourier_codes(V,np.random.default_rng(0)),2))

import sys, json
from pr import *
for f in sys.argv[1:]:
    d = json.load(open(f)); runs = d if isinstance(d, list) else [d]
    for run in runs:
        r = run['r']; PR = PRmat(r)
        S, ubS, stS = mis(lambda a, b: PR[a][b] >= 1)
        Q, ubQ, stQ = mis(lambda a, b: PR[a][b] == 2)
        print(f, 'alpha(Gamma_E)', len(S), ubS, stS, '| alpha(Gamma_<=1)', len(Q), ubQ, stQ, '=> C >=', len(S) + len(Q), flush=True)

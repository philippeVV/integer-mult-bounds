import sys, json, os, glob, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
for f in sys.argv[1:]:
    name = f.split('_')[1]
    G = json.load(open(os.path.join(HERE, '..', 'expXP', 'graphs.json')))[name]
    EV = np.array(G['Pts'] + G['Lns']); V, d = EV.shape; k = d // 2
    OR = json.load(open(os.path.join(HERE, f'orient_{name}.json'))); TAIL, HEAD = np.array(OR['tail']), np.array(OR['head'])
    Q = np.load(f); m = []
    for v in range(V):
        li = [i for i, e in enumerate(EV[v]) if HEAD[e] == v]; lo = [i for i, e in enumerate(EV[v]) if TAIL[e] == v]
        m.append(min(abs(np.linalg.det(Q[v][li])) ** 2, abs(np.linalg.det(Q[v][lo])) ** 2) * 2)
    m = np.log(np.array(m)); w = np.load(f.replace('sol_', 'solw_'))
    print(f, 'raw io min %.1f  #codes raw<e^-10: %d' % (m.min(), (m < -10).sum()), 'eig counts <1e-12,1e-9,1e-6:', [int((w < t).sum()) for t in (1e-12, 1e-9, 1e-6)])

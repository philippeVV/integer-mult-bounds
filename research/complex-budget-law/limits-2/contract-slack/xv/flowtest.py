import json, numpy as np, sys
from flowlp import fstar
name = sys.argv[1]; S = int(sys.argv[2])
G = json.load(open('../expXP/graphs.json'))[name]; EV = G['Pts'] + G['Lns']; E = G['E']
ends = [[] for _ in range(E)]
for v, es in enumerate(EV):
    for e in es: ends[e].append(v)
ends = np.array(ends); OR = json.load(open(f'orient_{name}.json'))
rng = np.random.default_rng(0)
for t in range(2):
    Rs = rng.choice(OR['R'], S, replace=False)
    print(json.dumps(fstar(ends, Rs, len(EV))))

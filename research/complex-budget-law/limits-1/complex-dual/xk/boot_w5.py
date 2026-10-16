import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'xg'))
from boot import boot
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from xk_sparse import load_any
for name in ['W3', 'W4M', 'W5']:
    V, nP, E, g = load_any(name); c, o = boot(V, E, R=60)
    print(json.dumps(dict(graph=name, E=E, deg=len(V[0]), boot=c, ratio=round(c / E, 5), threshold=1 / g)), flush=True)

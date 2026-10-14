import gzip, json, os
HERE = os.path.dirname(os.path.abspath(__file__))
def dot(a, b): return bin(a & b).count('1') & 1
def load():
    d = json.load(gzip.open(os.path.join(HERE, 'src/candidate/strict/global/base/providers/gauge44/CANONICAL-WORD.json.gz')))
    return d['ports'], d
P, W = load()
def rchoices(q): return [r for r in range(512) if dot(r, q) == 0 and dot(r, r) == 1]

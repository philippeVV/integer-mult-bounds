import gzip, hashlib, json, os, tarfile
HERE = os.path.dirname(os.path.abspath(__file__))
def dot(a, b): return bin(a & b).count('1') & 1
def load():
    # selected canonical word, read straight from main's pinned replay bundle
    tar = os.path.join(HERE, '..', '..', 'community-frontier1094-audit', 'selected-replay.tar.gz')
    with tarfile.open(tar) as t:
        raw = t.extractfile('candidate/strict/global/base/providers/gauge44/CANONICAL-WORD.json.gz').read()
    assert hashlib.sha256(raw).hexdigest().startswith('16bc0ebb04b4be6c'), 'unexpected canonical word'
    d = json.loads(gzip.decompress(raw))
    return d['ports'], d
P, W = load()
def rchoices(q): return [r for r in range(512) if dot(r, q) == 0 and dot(r, r) == 1]

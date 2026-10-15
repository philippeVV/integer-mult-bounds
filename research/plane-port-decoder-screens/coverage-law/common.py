"""Shared loaders for expSB (coverage-law scouting). Reads the pinned selected word only."""
import gzip, hashlib, json, os, sys, tarfile
from fractions import Fraction as Q
MAIN = os.environ.get('REPO', '/path/to/main/checkout')
CEIL = os.path.join(MAIN, 'research/complex-local-word-ceiling')
sys.path.insert(0, CEIL)
import ceiling as C  # pinned, unchanged

WORD = {1: 2412, 2: 418, 3: 432, 4: 142, 5: 82, 8: 8}
M, FAM, ST = 73, 896, 9


def dot(a, b):
    return bin(a & b).count('1') & 1


def load_ports():
    tar = os.path.join(MAIN, 'research/community-frontier1094-audit/selected-replay.tar.gz')
    with tarfile.open(tar) as t:
        raw = t.extractfile('candidate/strict/global/base/providers/gauge44/CANONICAL-WORD.json.gz').read()
    assert hashlib.sha256(raw).hexdigest().startswith('16bc0ebb04b4be6c'), 'unexpected canonical word'
    d = json.loads(gzip.decompress(raw))
    return d['ports']


def e8map():
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lower', 'e8map.json')
    return {int(k): v for k, v in json.load(open(p)).items()}


def rank_f2(rows):
    """rank over F2 of a list of int bitmasks."""
    basis = []
    for v in rows:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return len(basis)


def span(vecs):
    basis = []
    for v in vecs:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return basis


def pair_rank(A, B):
    """F2 rank of the dot pairing between span(A) and span(B)."""
    a, b = span(A), span(B)
    rows = []
    for x in a:
        r = 0
        for j, y in enumerate(b):
            r |= dot(x, y) << j
        rows.append(r)
    return rank_f2(rows)


def cap(local, centre=64, unpaid=1, grid=Q(1, 10 ** 9)):
    return C.ceiling(local, M, FAM, ST, centre, grid=grid, unpaid=unpaid)[1]

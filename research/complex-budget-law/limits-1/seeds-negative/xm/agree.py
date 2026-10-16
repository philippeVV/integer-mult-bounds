"""limits-1: the two split primes P0, P1 (fam.py 0 / fam.py 1) give the same dim W on every exact family row."""
import json, os
H = os.path.dirname(os.path.abspath(__file__))
a, b = (json.load(open(os.path.join(H, f'fam_P{i}.json'))) for i in (0, 1))
assert len(a) == len(b) and all(x['fam'] == y['fam'] and x['S'] == y['S'] for x, y in zip(a, b))
print(f'{len(a)} rows, dims agree mod both primes; max S = {max(x["S"] for x in a)}')
print('PASS_EXPXM_AGREE')

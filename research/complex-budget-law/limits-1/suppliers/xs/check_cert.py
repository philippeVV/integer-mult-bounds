"""Independent stdlib checker for trace_{inst}.json.
Checks: each listed path is a simple path of 5 INTERNAL edges from c_k to p_k of its session; x and y lie on the listed
paths; |pos gap| differs between the two sessions.  By Lemma T (tr of each rank-<=1 step = 1 on an exact session's
5-path, sign sigma_k = +-1 for orientation) this proves: the two sessions cannot both have rank(G_p - G_c - I) = 0 while
every internal segment on those three paths has rank <= 1 (any field of char != 2, in particular Q and Z[1/2]).
Usage: python3 -I check_cert.py inst.json trace_inst.json
"""
import json, sys
inst = json.load(open(sys.argv[1])); cert = json.load(open(sys.argv[2]))['cert']
E = set()
for a, b in inst['internal']: E.add((a, b)); E.add((b, a))
R = [tuple(s) for s in inst['sessions']]
def pos_on(k, path, v):
    c, p = R[k]
    assert path[0] == c and path[-1] == p and len(path) == 6 and len(set(path)) == 6, 'bad endpoints'
    assert all((path[i], path[i + 1]) in E for i in range(5)), 'non-internal edge'
    return path.index(v)
x, y = cert['x'], cert['y']
k1, k2 = cert['k1'], cert['k2']
assert [list(R[k1]), list(R[k2])] == [list(s) for s in cert['sessions']]
g1 = pos_on(k1, cert['paths_k1'][1], y) - pos_on(k1, cert['paths_k1'][0], x)
g2 = pos_on(k2, cert['paths_k2'][1], y) - pos_on(k2, cert['paths_k2'][0], x)
print(f'session {k1}: t(y)-t(x) = sigma*{g1};  session {k2}: t(y)-t(x) = sigma*{g2}')
assert abs(g1) != abs(g2)
print('CERT_OK: no sign choice reconciles them (difference', abs(abs(g1) - abs(g2)), ')')

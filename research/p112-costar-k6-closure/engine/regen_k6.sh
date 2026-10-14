#!/bin/sh
# Regenerate and independently replay the k=6 certificates (not committed: certk6.pkl is ~490 MB).
#   cliques.pkl          all 262,928 orthogonal cliques of P112 (cliques.py)
#   certk6_filter.pkl    initial L1 filter over 12,528 covering orbits x 56 capacity compositions
#   certk6.pkl           engine4 certificates for the 63,824 filter survivors
#   k6chunks/c00..c39    split6.py: 701,545 replay jobs; k6_open.pkl: the 23 open cases
#   replay6.py           verify_cert.replay on every job (independent of model.py / engine*.py)
# JOBS (default 4) sets the worker count.  Run under a MemoryMax scope (6 GB sufficed in a partial rerun).
set -eu
cd "$(dirname "$0")"
if [ -e k6chunks ]; then echo "k6chunks/ exists: remove it first (replay6.py skips chunks with a .ok marker)"; exit 1; fi
python3 cliques.py
python3 run_k6b.py | tee run_k6b.log
cp k6_open.pkl k6_open.committed.pkl
mkdir k6chunks
python3 split6.py
if cmp -s k6_open.pkl k6_open.committed.pkl; then echo "k6_open.pkl: identical to committed copy"; else echo "k6_open.pkl: DIFFERS from committed copy"; fi
mv k6_open.committed.pkl k6_open.pkl
python3 replay6.py
(cd k6chunks && sha256sum c??.pkl) > k6chunks.sha256
echo "sha256 of k6chunks.sha256: $(sha256sum k6chunks.sha256 | cut -d' ' -f1)"

#!/bin/sh
# Run from the repository root (Python 3.10+, stdlib only; ~1-2 min).
set -e
D=research/joint-query-first-pair-catalogue
AR=research/community-frontier1042-audit/selected-replay.tar.gz
python3 -I -B $D/pair_catalogue.py $AR research/joint-query/evidence/QQ16-SELECTED-WITNESS.json ${OUT:-/dev/null} | tail -2
python3 -I -B $D/trajectory_check.py $AR $D | tail -1
for c in ctl_47 ctl_noA; do
  if python3 -I -B $D/$c.py $AR $D >/dev/null 2>&1; then echo "CONTROL $c NOT REJECTED"; exit 1; else echo "control $c rejected"; fi
done

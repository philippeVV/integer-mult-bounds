#!/bin/bash
# Lemma S on all stage-2 survivors, --frames, resumable.
# usage: s3run.sh PARITY   (0|1; run two copies, one per parity, for 2 workers)
# Splits cases_k7_all.json into s3/kNN.json (100 cases each, in file order) if not yet done,
# then runs check_S7.py --frames on every s3/kNN.json with NN % 2 == PARITY.
# env PY: Python command with python-sat installed (default python3).
cd "$(dirname "$0")"
PY=${PY:-python3}
mkdir -p s3
[ -f s3/k00.json ] || $PY -c '
import json
a = json.load(open("cases_k7_all.json"))
for i in range(0, len(a), 100): json.dump(a[i:i + 100], open("s3/k%02d.json" % (i // 100), "w"))'
for f in s3/k*.json; do
  n=$(basename $f .json); n=${n#k}; [ $((10#$n % 2)) -eq $1 ] || continue
  [ -f s3/k$n.done ] && continue
  $PY check_S7.py $f --frames > s3/k$n.log 2>&1 && grep -q "LEMMA S:" s3/k$n.log && touch s3/k$n.done
done

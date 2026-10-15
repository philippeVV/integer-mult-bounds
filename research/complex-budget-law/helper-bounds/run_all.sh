#!/bin/bash
# Reproduce (about 1 minute). Needs PKG=<pr2474 checkout>/research/gmx-t74-borrow43-handoff and REPO=<main checkout>.
set -e
P=$PKG/new/borrow/words
${PY:-python3} -B convert_raw.py $REPO/research/k8-rank7-composition/sources/K8-SLOT-WORD.json.gz sel1094.json.gz
for x in base extra; do
  ${PY:-python3} -B liveness.py $P/BORROW-WORD-43$x.json.gz --json out/live-$x.json > /dev/null
  ${PY:-python3} -B chains.py   $P/BORROW-WORD-43$x.json.gz 2>/dev/null > out/chains-$x.json
  ${PY:-python3} -B propbound.py $P/BORROW-WORD-43$x.json.gz > out/prop-$x.json
  ${PY:-python3} -B propdag.py  $P/BORROW-WORD-43$x.json.gz 2>/dev/null > out/propdag-$x.json
done
${PY:-python3} -B liveness.py sel1094.json.gz --json out/live-sel.json > /dev/null
${PY:-python3} -B pinbound.py sel1094.json.gz --json out/pin-sel-rows.json > out/pin-sel.json
${PY:-python3} -B propbound.py sel1094.json.gz > out/prop-sel.json
CONE_ONLY=1 ${PY:-python3} -B propbound.py sel1094.json.gz > out/prop-sel-cone.json
${PY:-python3} -B propdag.py sel1094.json.gz 2>/dev/null > out/propdag-sel.json
echo done

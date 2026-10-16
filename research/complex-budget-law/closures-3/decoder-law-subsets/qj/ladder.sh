#!/bin/bash
# expQJ: rebuild the deletion-threshold table and rerun the low-root witness branch and bound (#968 ladder, p = 25..41).
# Build lines (binaries go to a temporary directory, removed afterwards):
#   gcc -O2 -o sub sub.c ; gcc -O2 -o cand cand.c ; gcc -O2 -o bb2 bb2.c
# The witness file (83 MB, 1,943,160 lines) is regenerated, not shipped:
#   python3 -I -B mkall.py allroots_in.txt && ./cand allroots_in.txt 41 | awk '$2==6||$2==7' | sort -u > w67.txt
# usage: bash ladder.sh [--full]   (default p = 25..36, about 1 min; --full adds p = 37..41, about 12 min more)
set -u
cd "$(dirname "$0")"
H=$PWD; B=$(mktemp -d); trap 'rm -rf "$B"' EXIT
for x in sub cand bb2; do gcc -O2 -o "$B/$x" $x.c || { echo FAIL_EXPQJ_LADDER build; exit 1; }; done
cd "$B"; ok=1
python3 -I -B "$H/mkin.py" sub_in.txt && ./sub sub_in.txt > sub_out.txt 2> /dev/null
python3 -I -B "$H/thr.py" sub_out.txt > thr.txt
if cmp -s thr.txt "$H/thresholds.txt"; then echo "deletion-threshold table: identical to thresholds.txt"; else echo "threshold table differs"; ok=0; fi
python3 -I -B "$H/mkall.py" allroots_in.txt && ./cand allroots_in.txt 41 | awk '$2==6||$2==7' | sort -u > w67.txt
echo "witness lines: $(wc -l < w67.txt) (expected 1943160)"; [ "$(wc -l < w67.txt)" = 1943160 ] || ok=0
PMAX=36; [ "${1:-}" = "--full" ] && PMAX=41
for p in $(seq 25 $PMAX); do
    lim=$(( (p - 8) / 2 / 2 )); r=$(./bb2 w67.txt $p $lim); echo "$r"
    [[ "$r" == *EXCEEDS* ]] && ok=0
done
[ $ok = 1 ] && echo PASS_EXPQJ_LADDER || echo FAIL_EXPQJ_LADDER

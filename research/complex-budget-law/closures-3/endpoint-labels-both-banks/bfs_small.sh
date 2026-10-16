#!/bin/bash
# expQW: build bfs.c and replay the two smallest exhaustive both-bank instances (V = 1, one helper).
# Build line: gcc -O2 -o bfs bfs.c   (the binary goes to a temporary directory and is removed afterwards)
# Each instance asks whether some schedule has loss L <= sum u - 1; the expected answer is "none".
# The larger instances of bfs1.out (up to 9.9e7 states, V = 2) use the same commands with log2cap 28; see README.
set -u
cd "$(dirname "$0")"
B=$(mktemp -d); trap 'rm -rf "$B"' EXIT
gcc -O2 -o "$B/bfs" bfs.c || { echo FAIL_EXPQW_BFS build; exit 1; }
ok=1
for inst in "3 3 1 1,2" "4 2 1 1,2"; do
    set -- $inst
    python3 -I -B mkinst.py "$1" "$2" "$3" "$4" "$B/inst.txt" > /dev/null
    r=$("$B/bfs" "$B/inst.txt" 1 22 2>/dev/null)
    echo "h=$1 GF($2) H=$3 ports=$4: $r"
    [[ "$r" == none* ]] || ok=0
done
[ $ok = 1 ] && echo PASS_EXPQW_BFS || echo FAIL_EXPQW_BFS

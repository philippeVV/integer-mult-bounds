#!/bin/bash
cd "$(dirname "$0")"
QQ16_CONTRACT=final QQ16_LIB=1 timeout 7200 ../expO/run modelcheck.py src/selected-replay.tar.gz 3 11 > modelcheck_m3_final.log 2>&1 &
QQ16_CONTRACT=finalf QQ16_LIB=1 timeout 7200 ../expO/run modelcheck.py src/selected-replay.tar.gz 3 11 > modelcheck_m3_finalf.log 2>&1 &
wait

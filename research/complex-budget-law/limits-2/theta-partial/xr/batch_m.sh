#!/bin/bash
cd "$(dirname "$0")"
lane() { for job in "$@"; do set -- $job; timeout 1200 ../expO/run m_adv.py W3 $1 $2 $3 30 2>&1 | grep -E "SUMMARY|base|ray re-solve"; done; }
lane "sic 11 0.05" "sic 12 0.3" "sic 13 1" "sic 14 2" "4 11 0.15" "4 12 0.6" "4 13 1.5" "bdry0.9 11 0.1" "bdry0.9 12 0.5" "bdry1 13 0.3" > batch_mA.log &
lane "2 21 0.1" "2 22 0.4" "2 23 1" "bdry0.97 21 0.2" "bdry0.97 22 1" "0.3 21 0.3" "0.3 22 1" "1 23 0.5" "sic 24 0.6" "4 24 3" > batch_mB.log &
wait

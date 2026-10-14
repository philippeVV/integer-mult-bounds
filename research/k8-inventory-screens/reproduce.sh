#!/bin/bash
# Reproduce both screens from a rebuilt k8 leaver chain.  Run from the repository root:
#   research/k8-inventory-screens/reproduce.sh PIPE_DIR PR635_DIR OUT_DIR
# PIPE_DIR must hold the outputs of README step 1 (notri-base, lv, lvc).  About 3 minutes after step 1.
set -euo pipefail
PIPE=$1; PR635=$2; OUT=$3; H=research/k8-inventory-screens; CC=research/gauge44-leaver8/construction/leaver/CODED-CELLS-k8.json
mkdir -p "$OUT"
FINAL=$PIPE/lvc/frames/FULL-WORD.json.gz
python3 -c "import gzip,hashlib,sys;h=hashlib.sha256(gzip.open(sys.argv[1]).read()).hexdigest();assert h.startswith('db2c7cfb'),h" "$FINAL"
# A. K-block cold ends: per-block pairing, pigeonhole premises on three words, counterfactual
python3 -B $H/kblock/diag1.py "$PIPE/lv/hpe/FULL-WORD.json.gz" "$PIPE/notri-base/hamming" > "$OUT/diag1-lv.txt"
python3 -B $H/kblock/pigeonhole.py "$PIPE/notri-base/hamming" $CC "$PIPE/lv/raw/FULL-WORD.json.gz" \
  "$PIPE/lvc/hpe/FULL-WORD.json.gz" "$PIPE/lv/hpe/FULL-WORD.json.gz" > "$OUT/pigeonhole.txt"
python3 -B $H/kblock/hypo.py "$PIPE/lv/hpe/FULL-WORD.json.gz" "$PIPE/notri-base/hamming" $CC > "$OUT/hypo.txt"
for x in noprio coordprio revcoordprio noexch exch0123 exchfirstcopy; do
  python3 -B $H/kblock/run_variant.py --main . --pr635 "$PR635" --word-dir "$PIPE/notri-base" --out "$OUT/v-$x" \
    --xparam $H/kblock/x-$x.json > "$OUT/v-$x.log" 2>&1
done
python3 -B $H/kblock/run_variant.py --main . --pr635 "$PR635" --word-dir "$PIPE/notri-base" --out "$OUT/v-helpersfirst" \
  --cgx $H/kblock/cgx_leaver_helpersfirst.py > "$OUT/v-helpersfirst.log" 2>&1
# B. E1->E2 promotions: exact backward ceilings, the 32 legal promotions, audited words
python3 -B $H/e1promo/analyze.py "$FINAL" "$OUT/final-analysis.json" "$OUT/cones-final.json" > "$OUT/analyze.log"
if [ -n "${FRAME_MODEL:-}" ]; then  # optional independent cross-check against PR1148's frame_model.py
  python3 -B $H/e1promo/xcheck_ceilings.py "$FINAL" "$OUT/cones-final.json" > "$OUT/xcheck.log"
fi
for name in single474 all32; do
  if [ $name = single474 ]; then SPEC="474:1"; else SPEC=$(cat $H/e1promo/evidence/PROMOS-32.txt); fi
  python3 -B $H/e1promo/make_e1_recipe.py "$FINAL" "$OUT/cones-final.json" "$OUT/RECIPE-E1-$name.json" $SPEC
  python3 -B research/gauge44-leaver8/construction/promo/build_promotion.py --input "$FINAL" \
    --recipe "$OUT/RECIPE-E1-$name.json" --out "$OUT/word-$name"
  python3 -B research/gauge44-leaver8/independent/audit_local.py --word "$OUT/word-$name/FULL-WORD.json.gz" \
    --main . --output "$OUT/AUDIT-$name.json" > "$OUT/audit-$name.log"
done
python3 -B $H/e1promo/neg104.py "$FINAL" "$OUT/cones-final.json" "$OUT/NEG-104.json" > "$OUT/neg104.log"
echo REPRODUCED "$OUT"

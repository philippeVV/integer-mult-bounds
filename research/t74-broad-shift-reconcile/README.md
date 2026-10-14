# T74 broad-shift reconciliation (review evidence for PR1897 / PR2044)

**This is an independent exact check. It discharges one review obligation and makes no κ claim.**

The maintainer review on PR1897 (`issuecomment-6174289936`) asked us to "reconcile the broad-shift guard correction with the inherited outer argument". The check is done here against `outer/GAUGEFREE-INTERIORS.md` and `check_gaugefree.py` at PR2044 head `322949a5`. The script uses sympy 1.13.3 and Fractions, and imports nothing from the PR.

## Findings

1. **The disabled `or True` assertion is outside the proof's branch.** `check_gaugefree.py:23` tests the shift s = −δ/λ. The proof uses s = −δ(1−θ), and the assertions for that shift are live.
   - Corrected, the inherited expression has the residual `−2bδ(2θy − δ(1+θ)²)/(1+θ)`.
   - This matches the correction in PR1876/PR1897.
2. **The inherited bound |δ| ≤ 1/2 holds, in fact strictly.** δ = β′_m/ρ, so |δ| ≤ 1/(2ρ). With the proof's shift, 1 − 2|s| ≥ 2θ/(1+θ) ≥ θ, which is what accretivity needs.
3. **GAUGEFREE-INTERIORS.md:32 is inaccurate.** It says |s| "can exceed 1/2 by θ/2" for s = −δ/λ. In fact |s| = |β′_m| ≤ 1/2, with equality at β′ = −1/2, where φ(π) = 0. So shrinking the shift is still needed, and the conclusion stands.
4. **The sentence "linear leftover has logarithmic range ≤ 2 b0 θ B (total)" is false.** It appears at :31 and in the analogous `2 b0 θ B'` at :84.
   - Exact counterexample: θ = 1/40, c = −19/12, B = 8, m = −48, δ = 203/492. The total range is 1421/3362·b0 ≈ 0.423·b0, which is more than 2θB·b0 = 0.4·b0.
   - The coefficient bounds stated at :31 only give ≤ 3 b0 θ B.
   - This is harmless: §5 (:92) charges 4 b0 θ B per diagonal, and the exact grid confirms that charge.
5. **Grid check.** It covers every reduced θ = p/q with q ≤ 40, c ∈ (1/12)ℤ ∩ [−2,2), every phase cell, every dyadic block and the centre anchor. That is 287,256 blocks and 7,568,034 exact entries of K = Drow·G_s·Dcol.
   - The §5 per-diagonal range bound θB²/4 + 4θB holds everywhere. The largest ratio is 1120/1681.

A side check enabled the second disabled assertion, `struct_check.py:130` (far blocks are zero), in a local copy. It passes exactly, so it hides nothing.

**Not addressed:** supplier admission, the outer schedule, the all-size, adapter and rounding obligations, and native reuse (O1, O2, O4–O6 in the reviews).

## Reproduce

```sh
uv run --no-project --python 3.12 --with sympy==1.13.3 \
  python research/t74-broad-shift-reconcile/broad_shift_reconcile.py 40 12
```

This takes about 3.5 minutes. It prints the JSON in `run_q40_c12.json` followed by `PASS_BROAD_SHIFT_RECONCILE`.

Credit: jeremykothe2 (the PR1673 → PR2044 line) and the PR1876 and PR1897 correction authors. Prepared by Philippe Van Velzen with substantial Claude (Anthropic) assistance. A sub-agent wrote the check, and the main session re-ran it with identical output. Apache-2.0.

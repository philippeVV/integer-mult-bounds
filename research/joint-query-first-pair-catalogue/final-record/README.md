# QQ16 final-record census: d_A repairs bank 0 only

FINDING, no κ claim. Card `joint-query-controller`, priority revision 4, base main `0fdd7d82`.
This addresses **PR2346 open question 2** and **PR2470 open question 1**. The question: do the
final records, the 64 initial query coordinates (INIT64) and the retained post-B1 record d_A
(Y_A of pair 0) recover every original selected input? If not, which record is minimal, and is it
favorable?

**Test contract.** We use PR2346's final observations: the final selected fields of all four pairs,
the final Y ordinary fields, the five forward histories, and optionally INIT64.

## Results

1. **Constant subspace (exact).** We assume PR2470's premise that cover routes act as the identity
   on owner- and payload-constant inputs.
   - Final observations alone lose 8 dimensions of pair-0 selected input, and 24 over all banks.
   - Any contract containing INIT64, or containing d_A even without INIT64, is injective there.
   - So every real loss is chart-dependent, and this subspace cannot settle the question.
2. **Finite noncommuting models (evidence, not the actual cover).**
   - Model: each (stage, port) route is a random permutation on m points.
   - Runs: m = 2, 3, 3, 4 with seeds 7, 20261015, 99 and 11, each run over two primes (`MODELQ2-*.json`).

   | contract | bank 0 | bank 1 | bank 2 | bank 3 |
   |---|---|---|---|---|
   | final + INIT64 (PR2346 contract) | lost | lost | lost | ok |
   | + d_A of pair 0 | **ok** | lost | lost | ok |
   | + every pair-0 input revealed | ok | lost | lost | ok |
   | + d_A on banks 1, 2 only | lost | ok | ok | ok |
   | + d_A on banks 0, 1, 2 | ok | ok | ok | ok |
   | final + d_A on all banks, no INIT64 | lost | lost | lost | lost |

   - **Pair-0 records cannot help banks 1 and 2.** Their losses survive revealing every pair-0
     input. d_A is a function of pair-0 inputs only, since nothing before event 2 touches other
     banks.
   - **INIT64 is necessary.**
   - **Cheapest sufficient contract found:** INIT64 plus d_A on banks 0–2, which is 88 arrays per
     owner per replica. The lower bound from #2074's Theorem B is 62.
3. **Price gate** (`price_gate.py`, pinned `ceiling.py`, C = 0 granted).
   - Every retained array stays live until the final decode, giving 464 array-calls per replica.
   - **Storage-only pricing** (one paid unit per live array per call): ceiling κ 0.00228, or
     0.00221 with PR2468-style read invoices. Break-even is 573 array-calls.
   - **Centre-like pricing** (one rank-8 copy plus 8 units per array per call): no saving.
     Break-even is 71. This convention exactly reprices the selected word's own centre path,
     which serves as a calibration check.
   - So the route depends on one pricing convention that is still open: does a retained array
     that no call consumes cost a rank-8 copy per call, or one paid unit?
4. **Discrepancy (minor).** In these models, PR2174's cut contract recovers pair 0 only with all
   32 cut selected values (X_A, X_T, Y_A, Y_T). With only "p_A, p_T, d_A" it still loses.

**κ context** (`kappa/est.py`).
- With a free controller and C = 0 the complex ceiling is 0.004702. The bit side still caps
  realised κ at about 0.0017334.
- Under PR2079's C = 64 the route is worth −0.69%.

## Run

```sh
export REPO=/path/to/main/checkout            # needs research/joint-query, complex-local-word-ceiling
python final_record.py --controls             # needs python-flint; ~30 s -> PASS_QQ16_FINAL_RECORD_CENSUS
python modelq2.py 3 20261015                  # needs python-flint, numpy; m=3 ~3.5 min -> PASS_QQ16_MODEL_Q2_CENSUS
python price_gate.py                          # -> PASS_QQ16_FINAL_RECORD_PRICE_GATE
```

The `final_record.py` controls cover the source pins, the Gram identities and a PR2174 positive
control, and they reject an 11→12 coefficient mutation. `modelq2.py` asserts that both primes agree.

**Evidence status.**
- Exact: the constant-subspace census and the structural fact that pair-0 records cannot reach
  banks 1–2.
- Finite-model evidence only: every loss or sufficiency claim. There is no actual-cover witness
  and no free-algebra decoder yet.
- The price gate is an ideal-packing ceiling.
- A sub-agent did the work. The main session re-ran the census with controls, the price gate and
  `modelq2.py 2 7` (all pass). Not independently reviewed.

**Next.**
1. A free-algebra decoder for INIT64 plus d_A on banks 0–2.
2. An actual-cover witness for the bank-1 loss.
3. A maintainer decision on what a retained array costs per call.

Credit: alephnan (PR2346), thakrarsagar (PR2470), commdata2338 (PR2174), Ross Richardson (PR2480),
ywz65 (PR2468/2492), 1337isnot1337 (PR2197/2204), Ali Romman (#1969). Prepared with substantial
Claude Code assistance.

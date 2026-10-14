# Two helper-inventory screens on the k8 leaver word

**Scoped negative results. No κ claim.** Both questions change the helper
endpoint inventory of the selected k8 leaver word (`research/gauge44-leaver8`,
final word `db2c7cfb…`), and both are answered with exact computations:

- **A. Cold ends in the leaver K-blocks.** No circuit on a leaver cell's own
  registers can restore rank-four cold ends in blocks {0,1,2,11,12,13}. The
  obstruction is a short pigeonhole argument on the native frame rule.
- **B. E1→E2 entrance promotions.** Only 32 of the 136 E1 entering helpers can
  be promoted, and every one of those promotions makes the paid profile worse.

Card `complex-new-helper-inventory` (chooser seed
`3c585d1f4e91703422bf417f1ab21f40`, catalogue revision 2), refined to these two
helper types. Base: main `ef5f39403610d231a0c27e8a6b4f08c10209bcd6`.

## A. Cold ends in the K-blocks

On the no-triangle word, `cold44x.py` pairs at rank four in all 14 blocks, giving 56
cold ends. After #1005's eight mutable-leaver cells, the chain pins pairing ranks
`[0,0,0,4,4,4,4,4,4,4,4,0,0,0]`, giving 32. Recomputed exactly on both cold44 input
words (`f45d193d…`, `598543fa…`):

| Blocks | \|I\| | dim L | dim K | rank(LK) |
|---|---:|---:|---:|---:|
| 3–10 | 8 | 4 | 4 | 4 |
| 0, 1, 2, 11, 12, 13 | 4 | **0** | 4 | 0 |

In each K-block, `I` is exactly the four retired complement registers (block 0:
225, 227, 229, 234). Their clean forms are rows of the invertible completion
matrix, so `P[I, sources]` has rank 4 and `L = 0`. In the no-triangle word the L
rows came from four face-event operands. Those operands are now leaver sources,
and they end at rank 6–7, outside the block frame. The 24 helpers that #1005
removed were packet copies, not cold-end carriers.

**Lemma (pigeonhole).**
- **Frame rule.** For a gate `t += q·s`, the gate frame contains both operands'
  frames, and afterwards both operands carry it.
- **Setting.** Consider a cell with leaver sources `x_0..x_5` and targets
  `y_0..y_5` that must deliver `Σ_{j≠i} C_ij x_j` to `y_i`, with every
  off-diagonal `C_ij ≠ 0`.
- **Premises.** The following were checked for all eight cells, on the raw word
  and on both cold44 input words:
  - **P1.** The start frame of `x_i` is not inside the port hyperplane of `y_i`.
  - **P2.** A source's start frame lies in another source's block frame only if
    the two sources share a block.
  - **P3.** No gate outside the cell touches a leaver source after the cell
    starts.
- **Claim.** In any circuit on the cell's 12 registers, no leaver source ends
  inside its block frame.

*Proof.*
1. Let `g_i` be the last cell gate on `y_i`, and let `v_i` be its other
   operand. Afterwards `y_i` and `v_i` share a frame that contains every
   `x_j`, `j ≠ i`.
2. `v_i` is not a target `y_m`. Otherwise `y_m` would contain `x_m`, which P1
   forbids, so `v_i` is a leaver source.
3. If `v_i = v_k` for `i ≠ k`, the later of `g_i` and `g_k` joins its target to
   a register that already holds that target's own source. P1 forbids this
   too.
4. So `i ↦ v_i` is a bijection onto the six sources. Every source ends holding
   five sources' frames from all three blocks, and by P2 none ends inside its
   block frame. ∎

The proof covers every relabelling, gate order, coefficient choice (including
cancellation) and circuit length. Consequently `dim L = 0` in every K-block for
every completion priority, source-copy exchange set and carrier placement.

Even if the four no-triangle-low sources were kept inside the block frame, each
K-block would pair only at rank 3, so the counterfactual is 18 cold ends, not 24.

**Finite cross-checks.** All of these give `(|I|, dim L, dim K, rank LK) = (4, 0, 4, 0)`
in every K-block:
- completion priorities PRIO6, none, coordinate rows and reversed coordinate
  rows;
- three source-copy exchange sets;
- helpers-first carrier placement;
- all 720 paired cell relabellings with order changes that keep the exact
  linear map (864 candidate (relabelling, source) cases, none kept low).

A CP-SAT gossip model (8 rounds, x0 kept low) is infeasible. The lemma
supersedes these checks.

## B. E1→E2 entrance promotions

**Legality.** In `audit_local.py` and `build_promotion.py`, a promotion moves only
start frames and A/B/T frame annotations; the target response is unchanged. Every
frame constraint is an upper bound:
- monotone chains;
- A/B boundary pins;
- target port-hyperplane ends;
- `helper_final`;
- low-end cones at the end of B.

So the minimal monotone completion of a raised start is legal exactly when the new
start lies in the helper's **backward ceiling**, the meet of everything downstream.
Legality is therefore independent across helpers.

**Enumeration.** These are exact GF(2) results and are identical on the promo
(`84d4da2e…`), hybrid (`c836d892…`) and final (`db2c7cfb…`) words. There are 192
entering helpers; 136 of them are E1.
- **104 have ceiling equal to their entrance**, so no E2 exists. These are 74
  width-8 and 30 width-7 helpers. Target port-hyperplane ends alone pin 88 of
  them. The other 16 also need the helper low ends or their cones.
- **32 have dimension-2 ceilings**, giving exactly one E2 each: 474–477,
  486–489, 522–529, 534–541, 550–553 and 562–565. None reaches E3.

The computed ceilings equal PR1148's independent `frame_model` upper bounds on
all 5,406 gates. All 208 forced E2 words for the pinned helpers are rejected by
`audit_local.py`.

**Price.** Each of the 32 promotions reshapes a rank-3 climb and makes two
B-partners gain a rank-1 child. Per promotion:
- Δprofile `{1:+1, 2:+5, 3:−4}`, with width 8→7.
- At the packing-aware frontier dual, ΔRC = **+0.0023566**. This is a loss of
  about 0.0066% κ per promotion.

The all-32 word, profile `{1:2444, 2:578, 3:290, 4:170, 5:68, 8:8}` and widths
`{4:56, 5:32, 7:62, 8:90, 9:240}`, passes `audit_local.py`.

#1240's certified packing ceiling falls from the base `0.001689769942883` to
`0.001683791850650` (−0.354%). The single promotion 474 gives `0.001689582481026`
(−0.011%).

An exact frame re-optimisation with PR1148's `frame_model` LP, integral here, gives
the same deltas: the promotions are additive and do not interact.

## Not covered

- **A:**
  - Gates across cells, or registers borrowed from outside a cell. None of the
    672 cell-target hyperplanes contains a block frame, but this is not a proof.
  - Changing non-cell factory parts.
  - Re-adding helpers, which undoes #1005.
  - Other frame rules, and rank-four frames other than block frames.
- **B:**
  - Changing the scalar program, target levels, centre frames, `helper_final`
    or end depths.
  - PR1148's own rank-seven-end word.
  - Global repacking, which was not run because both screens price negative.

Both screens are exact inside the native frame model for the fixed scalar word.
Neither is a theorem about other words, contracts or frame rules.

## Reproduce

1. Rebuild the k8 chain exactly as in `research/gauge44-leaver8/README.md`:
   - pinned PR635 `a6895a69…` → `reproduce.py --until word`;
   - `notri/reproduce_notri.py --word-only`;
   - `leaver/reproduce_leaver.py --id k8 --word-only`;
   - `leaver/reproduce_leaver_chain.py --id k8`.

   Put the outputs in `PIPE/{notri-base,notri-word,lv,lvc}`. This takes about 25
   minutes and reproduced every pin here, including the cell word `58e14a02…` and
   the final word `db2c7cfb…`.
2. From the repository root, run
   `research/k8-inventory-screens/reproduce.sh PIPE PR635_DIR OUT` (a few
   minutes). It writes fresh copies of the files under `kblock/evidence/` and
   `e1promo/evidence/`.

Optional extra runs:
- Set `FRAME_MODEL` to PR1148's standalone
  `research/k8-eight-ends/discovery/frame_model.py` (tomdif/integer-mult-bounds
  at `506ea08dd77388105682044ebb0848ce9b274964`) before running `reproduce.sh`.
  It cross-checks every gate ceiling; all 5,406 of 5,406 agree.
- `kblock/relabel_fast.py CODED-CELLS-k8.json` (about 12 minutes).
- `e1promo/run_reopt_all.sh FINAL-WORD OUT`. It needs Python 3.12 with SciPy via
  `uv`.

## Tool notes

- `promo/promo_controls.py` hard-codes E4/E5 dimensions, so on E1 recipes it
  reports failures instead of rejections. A generalised copy is in
  `e1promo/promo_controls_e1.py`.
- The `leaver/build_frames_lv.py` docstring names the k5 recipe.
- `cold44x.py` logs `len(K)` under the label "K". That value is `|I|`; dim K is 4
  in all 14 blocks.

## Open next questions

- **Cross-cell leaver circuits.** Could a gate path through another cell or a
  borrowed register keep one leaver source low?
  - Why it matters: each restored K-block direction is worth about +0.08% of
    the complex saving.
  - First step: extend P3 to allow the cell's sources to touch one borrowed
    register, and recheck the pigeonhole.
  - Blocked by: a borrowed final partner itself ends high, so the gain may be
    illusory.
- **Unpinning E1 ceilings.** 88 E1 ceilings are pinned by target
  port-hyperplane ends. Promotions therefore need changed target levels, which
  is a different contract (compare `coverage-law-change` in #1272).
  - First step: price a single target-level relaxation with the frontier dual.
  - Blocked by: the inherited data endpoint contract.

## Credit

- icekylinx (#645): gauge44 and the cold-space construction.
- PartialK: no-triangle word and PRIO6.
- jeremykothe2 (#1005): k8 chain with LeaverSAT/LeaverSAT2 cells.
- BlockInvert: `build_promotion.py`.
- tomdif (#1148): `frame_model` and the E4/E5 search.
- higherbeing (#1529): pairing-collapse observation.
- woklerx (#1126, #1138): leaver relabelling work.
- Douglas Colkitt (#1240): packing ceiling.
- Ali Romman (#1034): catalogue dual.

Prepared by Philippe Van Velzen with substantial Claude (Anthropic) assistance.
Parallel Claude Code sub-agents ran the experiments, and the main session
re-checked the audits, the profile deltas, the ceiling prices and the proof.
Apache-2.0.

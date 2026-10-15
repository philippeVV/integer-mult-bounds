# Helper-inventory lever: how far deletion can go, and proved frame-forcing bounds

Follow-up to `../README.md` and `../levers/`. This is a FINDING with no κ claim. It covers PR #2474's
borrowed words (head 5298f86c) and the selected #1094 k8 word.

**Model.** `wordmodel.py` re-implements the `audit_local` frame replay literally. As a control it
reproduces each word's `local_profile` and `width_inventory` exactly. On the selected word the
helper log-mass is 13802.03, which matches the budget law. The base and extra words give identical
numbers.

## (A) Alias headroom on #2474 (`chains.py`)

Helper streams can follow a retired register when their lifetimes do not overlap. The minimum number
of helper registers is streams − maximum matching − 112 source heads. Freed mass is priced as rank-48
idle with the package root (W = 66304, 666 calls), which reproduces κ = 0.0017036763879937.

| class | min helper registers | further deletions | priced headroom |
|---|---|---|---|
| fixed order, fixed gate frames (any lender, free rank-0 warps) | **422** of 429 | **≤ 7** | **≤ +0.381%** (κ ≤ 0.00171016), ceiling ×1.0040 |
| fixed order, frames re-solved (liveness only) | 336 relaxed / 344 literal | ≤ 93 / 85 | envelope +10.3% / +8.8% |
| reordering + frames re-solved | no liveness bound | up to all | ideal ×2.39 |

- **The fixed-frame rows are limited by the frame rule, not by liveness.** Every possible new alias
  is a helper whose first touch is at dimension 4–5, with a lender that retires at dimension 4.
  Exactly 9 such streams are unaliased (606, 635, 638, 639, 653, 667, 669, 670, 671), and at most 7
  can be matched.
- **The 240 width-9 helpers can never alias.** They start at dimension 0, and no register retires at
  dimension 0.
- **The reordering row is vacuous.** The maximum clique of streams forced to overlap under all
  reorderings is just the 112 sources.
- **Scope.** Only aliasing onto retired registers is counted; nested borrows inside gaps are not.

## (B) Lower bounds on helper frame-chain log-mass

**Lemma B1** (`propbound.py`): fixed gate list and order, any frame refit.
- **Bounds.** lo(g) is the span of the partners' static lower pins (source port lines, start pins,
  centre frames, target levels) together with lo of the previous gate on each register. hi(g) is the
  meet of the upper pins (target hyperplanes p_t^⊥, end pins, levels/centres) together with hi of the
  next gate. Warps act as barriers. Then lo(g) ⊆ F_g ⊆ hi(g) in every refit.
- **Proof.** Induction along each register's gate sequence: frames are monotone, and a gate
  equalises its registers' frames.
- **Cost bound.** A dimension DP with f(r) = r ln(73/r), which is concave so the DP is exact, gives
  each helper's minimum cost.
- **Control.** The actual frames lie inside [lo, hi] at every gate.

| | selected | #2474 (T74) |
|---|---|---|
| log-mass between first and last touch forced | 89.6% | 92.1% |
| full chains with start/end pins forced | 97.2% | 97.7% |
| ⇒ any frame refit gains at most | ceiling ×1.0188 | +1.74% (ceiling ×1.0174) |
| pure data-flow cone only (source lines forward, target hyperplanes backward) | 83.5% | 82.7% |

**Lemma B2** (`propdag.py`): any reordering that keeps non-commuting gates in order, plus any refit.
- **Method.** The same propagation runs over the dependency order. Each helper's span is then at
  least max dim lo − min dim hi.
- **Forced mass.** 44.9% (selected) and 45.6% (T74) of the first-to-last log-mass is forced.
- **Consequence.** The helper lever is capped at **×1.65** (selected; ideal ×2.46) and **×1.59**
  (T74; ideal ×2.38) unless the gate multiset changes.

**Where a map-level theorem breaks down** (`pinbound.py`).
- **Static pins alone force almost nothing.** Without coupling, they force 2.5% (selected): 444 of 480
  helpers are single-frame as far as their own partners' pins go. So "a target pin forces k climbs on
  a helper" is false helper by helper.
- **The 8D sub-lemma holds.** A helper read into target t at an 8D frame is stuck at p_t^⊥ and pays an
  8→9 climb. This happens for 210 streams on the selected word and 201 on T74.
- **What a map-level theorem would need.** The forcing comes from data flow through the gate DAG. A
  theorem about the map itself must bound the causal-cone ranks of every gate list that implements
  the decoder. That is not proved.

## Run

```sh
PKG=<pr2474>/research/gmx-t74-borrow43-handoff REPO=<main checkout> \
PY="python3" ./run_all.sh     # about 1 min. chains.py needs scipy; e.g.
# PY="uv run -q --no-project --python 3.12 --with numpy --with scipy python"
```

The main session re-ran everything from this directory, and every output in `out/` is byte-identical
to the sub-agent's run.

Prepared with substantial Claude Code assistance. A sub-agent did the work; it has not been
independently reviewed.

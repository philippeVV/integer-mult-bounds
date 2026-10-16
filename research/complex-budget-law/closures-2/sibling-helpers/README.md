# Sibling calls sharing helper registers: capped at ×1.41 (perfect joint code); realised ×1.00

Status: scoped cap with an open residual. Author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`.
Multipliers are over κ = 0.0017.

## Verdict

| Item | Over 0.0017 |
|---|---|
| Realisable today | **×1.00** |
| Ceiling, k = 2 sibling calls, perfect joint code on both suppliers | **×1.41** (×1.39 on T74; ×1.33 once one rank-1 junk buffer per shared helper is charged) |
| Ceiling, k = 4 | ×1.78 |
| Ceiling, helpers free (k → ∞) | ×2.43 |
| Complex side only | ×1.018 (the bit endpoint binds) |

The lever is capped, not closed. Readings (a) and (b) below are already used or capped. Reading (c) has no code, and no bound
excludes one. P(≥ +5% on both sides) ≈ 0.02 (judgement).

## Changed hypothesis tested

Today every call pays its own helper climbs:
- 3,730 of the 5,586 local mass on complex #1094;
- 3,454 of 5,310 on T74;
- 3,512 of 5,129 on bit #589.

The proposal: k sibling calls jointly use one helper pool, so each pays about 1/k of the helper children.

| Reading | Verdict |
|---|---|
| (a) time-sharing the registers | Already in today's word through the hosting lemma. The bill is coordinate advances, not storage. |
| (b) borrowing a sibling's dirty helper | This is restoration deferral across calls (`t37-restoration-deferral`, `z4-*`). It is mass-neutral, ≤ +2.81% complex, and ×1.018 overall. Siblings also sit in disjoint stage windows. |
| (c) joint use: one register carries forms of both calls | The only reading with a lever. It needs a merged same-stage call on two replicas **and** a vector-linear (rate-½) joint code on the decoder relay network. |

**Premise correction.** Helpers carry decoder data, not just dirt. In the realised network `absorb10`, 1,810 of the 2,553
register edges are helper steps. They deliver all 5,824 decoder pairs and the 448 centre pairs.

## Price (`price.py`)

**Model.** A fraction σ of the per-call helper children is served jointly by k calls. The freed coordinates go to ideal idle.

**Fitted law.** κ(σ,k) ≈ κ₀ / (1 − φ·σ·(1 − 1/k)), with:

| Word | φ |
|---|---|
| complex #1094 | 0.593 |
| complex T74 | 0.560 |
| bit #589 | 0.624 |

**Break-even, both sides, k = 2.**
- +5% needs σ ≥ 0.17;
- +10% needs σ ≥ 0.31;
- +25% needs σ ≥ 0.68.

## Lemmas (proved in text, not machine-checked)

**Lemma P (no mixing, no gain).** Suppose every shared register carries a form of only one copy at each step. Projecting onto
each copy then gives two valid single networks, so the joint cost is their sum. A gain needs genuinely mixed forms.

**Lemma S (reader cut).**
- Setup: a reader reads k shared registers holding forms (a_i, b_i), plus per-copy side information S_A ⊕ S_B.
- Statement: it gains at most k new pure-A plus pure-B dimensions in total.
- Proof: the A-pure and B-pure parts of W = span + S_A + S_B intersect in 0, so dim W_A + dim W_B ≤ dim W.
- Consequence: a shared-helper joint network needs |R(T′)| ≥ 2ρ(T′) for every target set T′. A single copy needs only ρ(T′).

**Non-discrimination of the old relaxation.** Under diagonal substitution (x_B = λx_A), every joint network yields a
single-copy network with 224 free relays. So any expSK-type lower bound for joint networks is ≤ LB(224 relays)/2 ≈ 70–140
per call. Such bounds can say nothing about this lever, so do not relaunch them.

## Exact census on `absorb10` (mod 2⁶¹−1)

`census.py`:
- every one of the 1,453 nodes has zero dependence slack (in-span rank 2,553 = in-edges 2,553);
- 112/112 targets need helper content;
- each target reads 3–8 helper registers, all with independent residuals.

`readercut.py`:

| Target set T′ | \|R(T′)\| − 2ρ(T′) |
|---|---|
| all 112 targets (\|R\| = 351, ρ = 50) | 251 |
| 90 shared-chain-node groups | min **0** (tight: e.g. 4 targets, 8 reads, ρ = 4) |
| greedy search, 16 seeds, ≤ 12 targets | min 1 |

No cut is violated, but several are tight. The skeleton has just enough reads for two copies.

## Stop test

Two stop tests were set. Neither has run or fired:
- **Prescribed (expSK MILP on two copies, stop if ≥ 95% of 1,607).** It cannot fire in this model class. The joint minimum
  is 137 per call, but that comes from the relaxation's non-discrimination, not from a code.
- **Decisive next test (not run).** Take one tight 4-target cone from `readercut.py` and ask exactly whether a shared code
  exists on the same helper edges, over F_p, solved bilinearly and then certified.
  - **Stop:** if no tight cone admits one, close the lever at σ ≈ 0.
  - **Continue:** if a cone works with ≤ 1.5× its single edges, price it at the measured σ and check the bit twin.

## Files and reproduction

Run each from this directory: `python3 -I -B <script>`. Everything is stdlib.

| Script | Runtime | Output |
|---|---|---|
| `price.py`: both pinned engines (`research/complex-local-word-ceiling`, `research/bit-revision3-levers-seed03`; root located relative to this file or by `IMB_REPO`) | 215 s | `price.out`, `PASS_EXPNU_PRICE` |
| `census.py`: per-node slack and per-target census | 0.2 s | `census.out`, `PASS_EXPNU_CENSUS` |
| `readercut.py`: Lemma S reader cuts | 79 s | `readercut.out`, `PASS_EXPNU_READERCUT` |
| `geo.py`: shared geometry. It reads `../../helper-resynthesis/expSK/{quads,alpha,instA}.json` and imports `../../helper-bounds/wordmodel.py`. | — | — |

`inputs/absorb10.json.gz` and `inputs/absorb10_net.json.gz` (48 KB) hold expSL's best realised helper network (R1+R2 design,
junk not charged).

The T74 row (`t74.py`, ×1.389 on the T74 root) and the two-copy MILP (`r1joint.py`, needs scipy and the expSK profiles) are
not included. Their outputs are quoted above from the scouting run.

## Scope

All prices are planning ceilings from pinned engines, not constructions. The census and reader cuts are exact on one realised
network, and they are necessary conditions only.

Prepared with Claude Code (Claude Opus 5.5) assistance. Not independently reviewed.

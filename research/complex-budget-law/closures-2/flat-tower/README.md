# Flat tower: more exchange pairs per replica (closed: realised ×1.00, ceiling ≤ ×1.05)

Status: scoped closure, author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`.
Multipliers are over κ = 0.0017.

## Verdict

| Item | Value |
|---|---|
| Realised | **×1.00** (κ unchanged). The selected 9-call / 4-pair skeleton is optimal among words with 2n+1 calls at m = 73. |
| Ceiling for anything left in this family | **≤ ×1.05** (κ ≤ 0.001789). It needs an unknown parallel word with e ≥ max(2, ⌈n/5⌉) extra units. Complex-only it is ×1.018, because the bit endpoint binds. |
| Correction to the earlier scouting row | The "variant F, n → ∞" row at **×1.32** (κ ≤ 0.002238) priced a flat tower without paying for its stages. It is **unrealisable**. Read it as ×1.00 realised and ≤ ×1.05 ceiling. |

## Changed hypothesis tested

The selected replica runs n = 4 exchange pairs with 2n+1 = 9 centre calls, and each call opens one fresh 8-block stage
(m = 1 + 8·9 = 73). The flat tower keeps m fixed and puts n > 4 pairs on one replica. To fit, two calls must share a stage.
There are two readings:

- **(a) same family:** a port register is routed through a stage it already used (the p→q→p return of the fivecall audit);
- **(b) different families:** two calls on different pairs sit at the same stage level, i.e. they run in parallel.

## Proof sketch

**Typed bound (`chain.py`).** Allow free block-diagonal bridges and paid L/U shears (rank r costs r; a mixed rank-1
transvection costs 2). Then n exchanges need ≥ 2n+1 centre units. The potentials a(·) = dim(·∩X) and b(·) = dim(·∩Y) must
each travel n, and the first pulse cannot be tight. The retained four-pair scalar chain attains 2n+1 for every n (checked
exactly for n = 1..8), and n = 4 is the selected nine-call skeleton (5 F, 4 B).

**(a) Same-family sharing needs a rank-8 repair (`stage.py`, exact GF(2) on the real ρ_(s,q) routes).**
- Fresh stage: all 896 (stage, port) pairs see the local view ⟨q⟩ and output W_s (+8).
- Used stage: all 4,032 (level, stage, port) routes and 1,440 p→q→p excursions see all of L9 (dimension 9), so no call is
  admissible.
- Route independence: for any g ∈ O(73,2), dim(gW ∩ L9) ≥ dim W − 64. A frame of dimension > 65 admits no call under any
  routing. The cheapest repair is a nested descent of rank 8, i.e. one full centre rank per port family.

**(b) Parallel calls are impossible at 2n+1 units (`levels.py`).**
- **Theorem.** Every typed word with exactly 2n+1 centre units has all blocks of rank 1. Let σ = n − a(S) − b(S) and
  τ = n − a(T) − b(T). With 2n+1 units, σ + τ ≤ 1 throughout, and every block rank is bounded by the slack (≤ 1).
- **Swap lemma.** Without descents, a call at a lower frame than the previous call's output commutes with it. The swap
  would make two same-direction shears adjacent and merge them into a rank-2 block, which the theorem forbids.
- **Consequence.** The frame rises by ≥ 8 per call, so 1 + 8·2n ≤ 65, which forces **n ≤ 4**.
- **With e extra units.** At depth 9 the word re-associates into ≤ 10 alternating blocks, so 2n + e ≤ 11e. This gives
  e ≥ max(2, ⌈n/5⌉) for n ≥ 5. It is a ceiling only, and no such word is known.

**Pricing (`price_flat.py`, `price_n.py`; pinned engines via `engines.py`).**
- (a) The proved repair mass is ≥ 326 per call at n = 5, growing to about 1,784. The coverage law allows a positive saving only
  below about 38–48 per call. At n = 5 the break-even repair is 2.16 mass per call, so the requirement exceeds it by ×151.
  The engines reject the word (a ≤ 0) even when the repair is priced as ideal idle.
- (b) With the counting bound on e, the ceiling is κ ≤ 0.001789 (×1.052), reached from n = 10.

## Stop test that fired

The test was: stop if a call that shares a stage needs a paid route repair of rank ≥ 8 (one centre). `stage.py` part B/C
proves exactly this. Then variant T applies, and `price_n.py` shows that n = 4 is the optimum of the tower
(n = 5, 6, 8, 16: ×0.97, ×0.95, ×0.90, ×0.73).

## Open residual

A depth-9 typed word with fewer than 9/4 units per pair (e ≥ 2) is not excluded. Its ceiling is ×1.05 (×1.018 complex-only).
Deciding it: is there a word at n = 5 with 12 units in ≤ 10 blocks? If not, retire the family outright. Low priority.

## Scope and caveats

- The bit side is assumed to have the same depth-9 structure (7-blocks, m = 64) and the same typed chain. The bit routes were
  not replayed.
- "Frame" means the binary C-frame of a family. Bridges are assumed to need equal frames, as the template's `add` enforces.
- Idle advances are treated as free beyond the coverage law, which favours the lever.

## Files and reproduction

Run each from this directory: `python3 -I -B <script>`. Everything is stdlib, and the repository root is located relative to
this file (override with `IMB_REPO`). The runtimes below are for Python 3.14, one worker, under a 4G memory cap.

| Script | Checks | Runtime | Output |
|---|---|---|---|
| `chain.py` | 2n+1 chain, n = 1..8, typed potentials | 0.1 s | `chain.out`, `PASS_EXPNL_CHAIN` |
| `levels.py` | σ/τ theorem inequalities, non-commuting calls, level assignment | 4.6 s | `levels.out`, `PASS_EXPNV_LEVELS` |
| `stage.py` | (a) on the 112 archived ports (reads `research/direct-delivery504-global/.../DATA-PHASES.json`) | 5.2 s | `stage.out`, `PASS_EXPNV_STAGE` |
| `price_flat.py` | (a) and (b) priced on both suppliers | 79 s | `price_flat.out`, `PASS_EXPNV_PRICE` |
| `price_n.py` | variants F and T for n = 2..1024, including the withdrawn ×1.32 row | 23 s | `price_n.out`, `PASS_EXPNL_PRICE` |
| `engines.py` | imports `research/complex-local-word-ceiling/ceiling.py` and `research/bit-revision3-levers-seed03/levers.py` unchanged | — | — |

The outputs match the original scouting logs byte for byte. A numeric corroboration (`merged.py`: n = 3 merged patterns stall
at residual 1.0 over 150 restarts, about 10 min) is not included. It adds nothing to the exact checks.

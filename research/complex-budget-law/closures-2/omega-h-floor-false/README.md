# The ω·h decoder-cost floor is false (exact counterexamples)

Status: a scoped refutation, with exact certificates. Author-checked, LOCAL-ONLY until published. No κ claim.
Base: main `0fdd7d82`. Multipliers are over κ = 0.0017.

## Verdict

**The conjectured floor "C ≥ ω·h for every masked decoder" is false in every form tried.** Here:
- ω is the clique number of the required graph {B_ab = 1};
- h = rank_F2 B;
- cost is Σ_k rank_F2 B[R_k, A_k] (the #446 Thm C cost model).

The floor was proposed as a general law and as a one-second family screen ("×10 at θ = ¼ needs V/ω ≥ 13"). **Neither use is
valid.** The counterexamples are punctures of E8/E7, so they are not a new compression mechanism, and the multiplier is
**×1.00**.

| System | V | ω | h | ω·h | C ≤ | |
|---|---|---|---|---|---|---|
| P120 (E8 lines) | 120 | 8 | 9 | 72 | 72 | equality |
| P119 = P120 − 1 line | 119 | 8 | 9 | 72 | **71** | violates |
| **P112** = P120 − an orthogonal frame (vertex-transitive) | 112 | 8 | 9 | 72 | **64** | violates |
| anneal subsets (`best_1.json`, `best_2.json`) | 74 / 73 | 8 | 9 | 72 | **56** | violates, with V ≥ ω·h |
| E7 (63 lines) | 63 | 7 | 7 | 49 | 49 | equality |
| E7 − 1 line | 62 | 7 | 7 | 49 | **48** | violates |
| **F_7** = all 64 odd vectors of F2^7 = E7 + all-ones | 64 | 8 | 7 | 56 | **50** | violates (closed system) |

Every decoder in the table is the frame-gauge G/2 restricted to the port set. It is checked entry by entry in exact
arithmetic, and the supports used in the cost are checked to be the true supports.

**Forms refuted.**
- literal (C ≥ ω·h);
- min form (C ≥ min(V, ω·h)), by P119, P112 and the anneal sets;
- closed-system form, by F_7;
- vertex-transitive form, by P112;
- the proposed denominator-2 absorption inequality: on P112 every G/2 term has denominator exactly 2, so N_int = 0 and the
  sum is 64 < 72.

**Correction to earlier notes.** P112 was listed as "attaining" ω·h. In fact its h is 9, not 8.

## Why it had to fail (proof sketch)

C is monotone under restriction, and deleting a port lowers C by at most 1. ω and h are locally constant under puncturing.
So any floor f(ω, h) valid for all systems is at most the cost of the smallest system with that (ω, h): a clique plus h − 1
rank-raising ports, so **f ≤ ω + h − 1**. The ω·h equalities on E8, E7 and BW16 are properties of those objects, not a law.

## What survives (proved, elementary)

**Monotone gain.** For Q ⊆ F, C(Q) ≤ C(F) ≤ C(Q) + |F∖Q|, so gain(Q) ≤ gain(F).

**Consequence.** A per-object proved floor C(F) ≥ c_F gives the family screen **γ(Q) ≥ 1 − (V_F − c_F)/|Q|** for every
Q ⊆ F. For example, C(P120) ≥ 72 would give γ(Q) ≥ 1 − 48/|Q|, which P104 attains. That is the only legitimate way to turn an
object floor into a family screen.

## Consequences for the map

- The "V/ω ≥ 13 for ×10 at θ = ¼" screen must not be used to close families. On record objects it is already wrong by
  ×1.13–1.29 (P112: screen .643 against .571 realised).
- The reading "BW16 has V = ω·h, so gain 0" (expL3) is void. BW16's window stays as stated, open and ≤ ×2.8.
- None of the map's proven large-port floors (γ₂ ≥ 0.161, γ₃ ≥ 0.357) relies on ω·h. They come from minrank/α bounds (L1′).

## Stop test

None was needed; this is a refutation. An exhaustive search at h ≤ 5 was requested but not run, because it is moot: the
floor already fails at h = 7 and h = 9. The open side question is the smallest window with positive gain (E7 at h = 7 is the
smallest known; h = 6 is open). It is not a κ lever.

## Files and reproduction

| Script | Run | Runtime | Output |
|---|---|---|---|
| `verify.py`: exact certificates for every row (stdlib) | `python3 -I -B verify.py` | 2.9 s | `verify.out`, `PASS_EXPQP_VERIFY` |
| `e8.py`: E8 line model and gauge decoders | — | — | — |
| `punct.py`: ω (required-clique search) and rank_F2 B | — | — | — |

Inputs:
- `best_1.json`, `best_2.json`: the anneal port sets and gauges.
- `e8map.json` (4 KB): a copy of `research/plane-port-decoder-screens/lower/e8map.json` from branch
  `findings/plane-port-decoder` (commit 9f245f32). It is not on main.

Not included: `p112t.py` (P112 vertex-transitivity: 40 random frame-stabiliser elements, one orbit of size 112), `e7f7.py`,
`closed.py` and `anneal.py`.

Prepared with Claude Code (Claude Opus 5.5) assistance. Not independently reviewed.

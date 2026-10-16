# Exact instances (expXQ, expXU): odd-q Hermitian schedules, and the frame screens rerun on the exact q = 9 instance

Status: exact constructions, exact checks and screens, author-checked, LOCAL-ONLY until published (the q = 9 schedule itself is
already public, see below). No κ claim. Base: main `0fdd7d82`. Multipliers are over κ = 0.0017. Inputs: `../../limits-1/suppliers/`
(expXJ: the odd-q Hermitian code and S(q); expXN: frame families), `../../limits-1/bootstrap/` (Boot, flow rule).

## Verdict

| Item | Value |
|---|---|
| Question (expXQ) | is the odd-q Hermitian Tanner code on the full PG(2,q), q = p² (S(q) = [(q+1)(q²+q+1) − √q(q²+1)]/4: 166, 3449, 26434), schedulable without helpers, i.e. is there an order with F(σ) = S? |
| **q = 9: EXACT** | the Singer cyclic order m = 16, c = 3 (P_i at 2(16i mod 91), L_l at 2((16l − 3) mod 91) + 1) has Boot cost 166 = S and a feasible back-arc transportation; the Eulerian 5-in/5-out orientation has \|R\| = 166, Γ − R acyclic; gate-by-gate over F₉ all 182 [10,5] gates are bijective and **rank(Φ − I) = 0**, so dim W = 166: exact and sound, no helper ports, no dropped sessions. This is Braverman–He's 166-session code (BH's own order needs 9 C + 9 J ports). Negative controls: perturbed line codes give rank 1 / 61, matching dense dim W 165 / 105 |
| Published | as the first exact helper-free schedule of the BH code, in a comment on [discussion #2838](https://github.com/CrocSwap/integer-mult-bounds/discussions/2838); artifact [`bh-exact-schedule` @ a5693d90](https://github.com/philippeVV/integer-mult-bounds/tree/a5693d90/research/coding-gap-exponent/bh-exact-schedule). The copies here are that artifact's `schedule9.json`, `check_json.py`, `verify.py` (with PASS markers and no rewrite, below) |
| **q = 25: excess +31** | best schedulable F = **3480 = S + 31** (cyclic, m = 13, c = 148); simulated over F₂₅: rank(Φ − I) = 31, dim W = 3449, an independent confirmation of S(25). Best Boot upper bound 3465 = S + 16, but every low-cost order found is not schedulable (true F 3528–4321) |
| **q = 49: excess +157** | best schedulable F = **26591 = S + 157** (cyclic, m = 185); Boot ≤ 26515 = S + 81 (Baer layout, own F ≥ 26758) |
| Trend | schedulable excess 0, 31, 157 ≈ 0, 0.99, 1.83 × q^{3/2}/4 (0.050 q² at 25, 0.065 q² at 49): growing in absolute terms, shrinking relative to S (0.9% → 0.6%). Boot excess 0, 16, 81 stays below q^{3/2}/4 (6.75, 31.25, 85.75). The cyclic order is a strict local optimum of F at q = 25 (118,590 exact F-evaluations around it, no improving move); orders with lower Boot cost break the transportation |
| **Frame screens on the exact q = 9 instance (expXU)** | all null: **345/345** GRS-F₂ ranks full (expXN's 23 phase-frame families × m = 4..8 × 3 char-2 draws); two-class theorem on **41/41** partitions (expXN's 35 + 6 cut-aligned), so nullity 0 for all X, Y and all m ≤ 15; **184** rational-frame rays (m = 2..8, 34 family directions from expQR and expXN), every charpoly gcd trivial and t = 1 nullity 0 |
| **Contract mismatch** | the schedule is a circuit over **F₉**, and F₉ has additive exponent 3: there is no injective additive map F₉ → F₂ᵏ, so no binary (XOR) encoding of its gates exists. Upstream `finite_bit_contract.py` accepts the incidence (166 roles, 182 gates, 1076 physical edges = 166 + 744 + 166) and **rejects** the circuit at the scalar check (`Wrong input-to-output permutation`). Not routable either: every head→tail path in Γ − R has 5 arcs and 5·166 = 830 > 744 forward arcs |
| Multiplier | **×0** realised; frame route **×1.0** on the exact instance (κ stays 0.0017). Conditional dimension value at q = 9 (expXJ): a(b = 5) = 0.068, κ ≤ 0.064 (×37), needing frames that the screens exclude |

## Reading

- The char-2 schedule obstruction (PG(2,16)−M: Boot ≤ 784 vs dim 772, `../../limits-1/bootstrap/`) **does not occur at q = 9 in
  odd characteristic**: the arithmetic-side Hermitian code at q = 9 meets expXJ's "schedule + soundness" condition exactly.
- At q = 25 and 49 no exact schedule was found; Boot = S is neither proved nor refuted there (the only lower bound is S itself).
  A dropped-session version exists (the cyclic schedule with 31 dead sessions at q = 25), a 0.9% / 0.6% session loss.
- The exact instance removes expXN's last caveat ("BH scenario A is not a full contract instance", the C/J lanes) and the frame
  route stays closed: the dense, full-rank char-2 and char-0 realisations of this incidence see no forced structure. Not covered:
  GL(m, F₂) address-map frames, and a family-free (universal) certificate.
- The contract admits this incidence only through a char-2 realisation (F₂ gates on the same roles and ports). Whether one exists
  is open (as for BH, expXN SAT undecided); if it existed, GRS-F₂ would already give s ≥ Wm for every phase frame on it.

## Next steps and stop tests

1. q = 25: an F-aware search that keeps the transportation feasible (reverse directed 4-cycles of the orientation, re-sort),
   seeded with the cyclic schedule. **Stop** if F stays ≥ 3480 after about 1 h.
2. Frame route on the exact instance: **none recommended.** If reopened: an F₂ realisation of the exact incidence
   (CryptoMiniSat with native XOR, W = 166, 4550 gate bits). Either answer leaves κ unchanged. Stop after 1 h undecided.

## Files and reproduction

| Script (run from the directory shown) | Checks | Runtime | Output |
|---|---|---|---|
| `xq/`: `python3 -I -B check_json.py schedule9.json` (stdlib) | perfect difference set, the 910 Singer incidences, Eulerian 5/5, Kahn acyclicity of Γ − R (without the stored order), \|R\| = S(9) = 166 | < 1 s | `check_json.out`, `PASS_EXPXQ_JSON` |
| `xq/`: `uv run -q --no-project --python 3.12 --with numpy python -I -B verify.py 3 cyc9.txt dense` | transportation 166/166, 182 gates bijective, rank(Φ − I) = 0, dense dim W = 166, `schedule9.json` byte-identical | 1 s | `verify9.out`, `PASS_EXPXQ_VERIFY9` |
| `xq/`: `PERTURB=one uv run … verify.py 3 cyc9.txt dense` | negative control: rank(Φ − I) = 1, dense dim W = 165 | 1 s | `verify9_one.out`, `PASS_EXPXQ_VERIFY9_ONE` |
| `xq/`: `uv run … mkgraph.py 3`, `… cyc.py 9` | graph files; Singer cyclic scan over all (m, c): best cost 166 = S | < 1 s | `PASS_EXPXQ_CYC9` |
| `xq/`: `uv run … mkgraph.py 5`, `… cyc.py 25` | best cyclic cost 3480 = S + 31 (class-sign orientation is Eulerian, so F = cost) | 2 s | `cyc25.out`, `PASS_EXPXQ_CYC25` |
| `xq/`: `uv run … verify.py 5 cyc25.txt` | the q = 25 cyclic schedule over F₂₅: transportation 3480/3480, 1302 gates bijective, rank(Φ − I) = 31 ⇒ dim W = 3449 = S(25) | 12 s | `verify25.out`, `PASS_EXPXQ_VERIFY25` |
| `xq/`: `uv run … mkgraph.py 7`, `… cyc.py 49` (slow) | best cyclic cost 26591 = S + 157 | 100 s | `cyc49.out`, `PASS_EXPXQ_CYC49` |
| `xu/`: `python3 -I -B summarize.py` (stdlib) | re-reads the logged expXU screens: 345 full GRS-F₂ ranks, 41 two-class partitions without subfield roots, 184 rays with trivial gcds | < 1 s | `summarize.out`, `PASS_EXPXU_SUMMARY` |

Changes against the #2838 copies: `verify.py` and `mkgraph.py` insert their own directory into `sys.path` (so `-I` finds the
local `code.py` and `field.py`, not the stdlib `code` module), `verify.py` compares with the stored `schedule9.json` instead of
rewriting it and prints a PASS marker (expected ranks 0 / 31 at q = 9 / 25; controls 1 and 61); `check_json.py` asserts and
prints a PASS marker; `cyc.py` compares with the stored `cyc{q}.txt` instead of rewriting it. `run_all.sh` deletes the
regenerated `graph49.txt` (1.1 MB). The q = 49 schedule is not simulated (26591 sessions); its excess is the exact F of the
cyclic order, and S(49) = 26434 is expXJ's exact value (`../../limits-1/suppliers/`, `odd_singer.py 7`).

Included but not replayed: `xq/`: `baer.py`, `quot.py`, `saF.py`, `classa.py`, `osa2.c`, `osa3.c`, `mcf.py`, `check_order.py`
(Baer layouts, quotient scans, SA on Boot cost and on exact F) with logs `lsa25_*.txt`, `lsaF25_1.txt`, `lcla25_F.txt`,
`lq*.txt`, `lver25_cyc.txt`, `baer9.txt`. `xu/`: `circ9.py` (instance; F₉, GF(2¹⁶) and mod-p transfers), `contract.py`
(upstream checker run; takes the repo path, default: the repository root found from the script path) with `contract.log`,
`contract9.json`; `screen9.py`, `twoclass9.py`, `twoclass_cut.py`, `ray9.py` with their logs and JSON results. They need
python-flint, and `ray9.py` imports expQR's `families.py`, which is not copied.

Sources: expXQ, expXU. Prepared with Claude Code (Claude Opus 5.5) assistance. The q = 9 schedule, the q = 25 rank and the
cyclic excesses are exact; the q = 25/49 searches are heuristic; the screens are exact ranks and gcds over the listed families.
Not independently reviewed.

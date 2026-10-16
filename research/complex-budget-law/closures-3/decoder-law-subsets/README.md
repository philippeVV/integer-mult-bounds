# Decoder law on subsets: 8-term decoders have γ ≥ ½ on every |Q| ≥ 79; rank-9 lifts closed on supports ≤ 26 lines

Status: scoped partial closure, author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`.
Multipliers are over κ = 0.0017. **Cross-posted to the decoder thread #1587**: this lane extends PR #968 (the 8-term
deletion ladder, `research/biframe-decoder-ladder`).

## Verdict

| Item | Value |
|---|---|
| Question (expNP Theorem W stop test) | is γ_min = min_{Q ⊆ P120} C(Q)/\|Q\| ≥ ½? The one-block column opens only if some Q has gain(Q) > \|Q\|/2 |
| 8-term (bi-frame) decoders | **γ ≥ ½ proved for every Q with \|Q\| ≥ 79** (p = \|P120 ∖ Q\| ≤ 41), extending #968's p ≤ 24. The independent checker covers p ≤ 36; p = 37–41 rest on the C code only |
| 8-term decoders, p ≥ 42 | **open**. Span-only relaxations provably fail for \|Q\| ≤ 47, and the decoder-free profile relaxation already fails at p ≥ 28 |
| Rank-9 lifts on non-clique / dependent supports (expL1 open case R3) | **Theorem R3**: under R2 and the spanning condition SC, the lift is a κ-family member, so credit ≤ 9 and cost ≥ 72 in every GL9 gauge. **SC holds for every support with ≤ 26 lines.** Residue R3″ (≥ 10 bad lines, smallest found \|S\| = 52) is open |
| Any decoder, \|Q\| ≤ 95 | **open**: needs a gain-decay statement gain(Q) ≤ \|Q\|/2 < 48 that the law C(P120) ≥ 72 does not imply |
| Cheapest known ratio | **P104, 56/104 ≈ 0.538** (#717/#968). Anneals find nothing below 0.543; at fixed sizes 96/88/72 the cheapest found are 57/57/56 |
| Multiplier | **×0** for the one-block column at every known decoder. No counterexample found anywhere |

## Tools (expNZ; proved, elementary)

Model: N = Σ_k p_k q_kᵀ with N_aa = 1 and N_ab = 0 on orthogonal pairs of Q; ρ_k = rank_F2 B[R_k, S_k], B = (J + G) mod 2.

- **(L1) Odd-join additivity.** If every pair across Q = Q₁ ⊔ Q₂ is odd, C(Q) = C(Q₁) + C(Q₂) (restrict each term; B[R, S] is
  block diagonal). So C is additive over components of the orthogonality graph.
- **(L2) Cluster floor.** A pairwise-orthogonal A has C(A) = |A|, so C(Q) ≥ cl(Q), the largest induced cluster subgraph. But
  **cl(P120) = 14** (CP-SAT optimal, `nz/cluster.py`), so this floor is useless for large Q.
- **(L3) Packing theorem (Cauchy–Binet).** The generic rank of Σ_k D(x_k) B′ D(y_k) equals the best packing of disjoint
  nonsingular blocks inside the supports. Consequences: every Lemma-B/L8 certificate is a combinatorial packing; and
  C(Q) ≥ |Q|/2 for **every** Q would follow from one rational B′ (support in diag ∪ orth, nonzero diagonal) whose
  Q-nonsingular blocks all have F2 rank ≥ half their size. Existence of such a B′ is open (0/1 and the zonal choices fail).
  The F2 certificate is tight on the P120 frame decoder: hill-climbing reaches 72 = cost (`nz/certhc.py`).
- **No finite certificate by enumeration:** ≥ 2^120/|W(E8)| ≈ 1.9·10^27 orbits; any proof must be algebraic.

## 8-term deletion ladder past p = 24 (expQJ)

**Threshold table** (`qj/sub.c`, every subspace of F2⁹, all 255 W-classes of supports, every span level; exact). Minimum
deletions for span(N(x) ∩ Q) ≤ s:

| support class | \|N\| | ≤ 8 | ≤ 7 | ≤ 6 | ≤ 5 | ≤ 4 | ≤ 3 | ≤ 2 | ≤ 1 |
|---|---|---|---|---|---|---|---|---|---|
| root, N(a) (E7) | 57 | 1 | **25** | 37 | 45 | 49 | 53 | 55 | 56 |
| D7 | 78 | 14 | **42** | 58 | 66 | 70 | 74 | 76 | 77 |
| next class | 83 | 27 | 51 | 63 | 71 | 75 | 79 | 81 | 82 |
| all other classes | ≥ 84 | ≥ 26 | ≥ 51 | ≥ 63 | ≥ 71 | ≥ 75 | ≥ 79 | ≥ 81 | ≥ 82 |

The ≤ 8 and ≤ 7 columns reproduce #968 exactly; ≤ 6 … ≤ 1 are new.

**Proof for 25 ≤ p ≤ 41.** With δ(x) = max(0, 8 − span(N(x) ∩ Q)), Lemma 2 gives β_k ≥ 7 − δ(v_k) − δ(w_k), so
cst ≥ 56 − Σ_k(δ(v_k) + δ(w_k)); |Q| ≥ 79 spans, so the decoder is a bi-frame.
1. For p ≤ 41 only root-type sides can be low (D7 needs 42 deletions, all others 51), and δ ≤ 2.
2. A root line is at most one v-side and one w-side, so Σ ≤ 2Δ(Q), Δ = Σ over root lines of δ.
3. Exact branch and bound over witness unions (`qj/bb2.c`, line 0 WLOG low by W(E8)-transitivity): max Δ = 1 (p = 25), 2
   (26–35), 6 (36–41).
4. So cst ≥ 54, 52, 44, which beats ⌈(120 − p)/2⌉ = 48, ≤ 47, ≤ 42. ∎

| p | proved floor | needed ⌈(120 − p)/2⌉ |
|---|---|---|
| 25 | 54 | 48 |
| 26–35 | 52 | ≤ 47 |
| 36–41 | 44 | ≤ 42 |

**Not closed: p ≥ 42.** D7 sides can lose span from p = 42 and all classes from p = 51; the witness count for one root grows
from 63 (p = 25) to 16,193 (p = 41). Any argument using only side spans plus thresholds is provably too weak: the decoder-free
profile relaxation fails from p ≥ 28 and has floor 0 for every |Q| ≤ 47. Small Q need true ranks and exactness.

**Numerics (not replayed):** the minimum cst found by annealing exact bi-frames is 56 from p = 16 to p = 68 (best margin 4 at
p = 15–17, ≥ 8 for every p ≥ 24); weight-only bases bottom out at 2c − n ≥ 20; frame decoders never go below 64.

**Next step.** Couple v_k with its dual w_k (a cheap term needs *both* sides low) and run the witness B&B on terms for
p = 42…70. *Stop:* a term-level witness union ≤ p with cost bound below ⌈(120 − p)/2⌉; realise it in Fractions.

## Theorem R3 (expR3): rank-9 lifts on dependent / non-clique supports

**Statement.** Assume R2 (X_T = [R_T, 0] on T = V ∖ S, rank 8) and SC (for every line a, the lines of orth(a) ∖ S span r_a^⊥).
Then the c_s are concurrent unless S is an A2 triple (closed in expL1), so the lift is a κ-family member with A = S, with
credit ≤ 9 and cost ≥ 72 in every GL9 gauge.

**Proof sketch.** At a good line a the condition x_s·y_a = 0 reads β_a c_s·r_a + η_a = 0, so (c_s − c_s′)·r_a = 0 for all
s, s′ ∈ S ∩ orth(a), at every good line, in T or in S. The lines orthogonal to {s, s′} span {r_s, r_s′}^⊥ (E6 or D6, rank 6), so
the affine lines c_s + K r_s meet pairwise; pairwise-meeting lines with distinct directions are concurrent or coplanar, and
coplanar with |S| ≥ 3 means the A2 triple. **SC threshold:** the rank-6 flats of E7 have at most 36 lines, so the minimum
cocircuit is 27, and SC holds whenever max_a |S ∩ orth(a)| ≤ 26, in particular for every |S| ≤ 26.

**Extension.** If SC fails at ≤ 9 lines, concurrency still holds and the Y side is an explicit lift; this closes all 7
minimal SC-violating (E7-cocircuit) types.

**Residue R3″ (open).** Supports with ≥ 10 bad lines; the smallest found has |S| = 52 (T = E7 ∪ {t0}, where R2 itself is
doubtful). Evidence against new components: at κ-points on these supports the moduli tangent dimension equals the κ-family
dimension (4/4).

**Checks.** Every exact member tested has gauge floor `mincost` 72 or 73; label enumeration on the true model gives credit ≤ 9;
every negative control exceeds its bound (dropping the forced pairings gives credit 10–12; dropping a D6 cocircuit gives extra
freedom; non-concurrent perturbations are infeasible).

## Files and reproduction

Runtimes under a 4G cap, one worker.

| Script (run from its subdirectory) | Checks | Runtime | Output |
|---|---|---|---|
| `qj/`: `bash ladder.sh` | builds `sub.c`, `cand.c`, `bb2.c` (`gcc -O2 -o sub sub.c` etc., in a temporary directory); rebuilds the threshold table (identical to `thresholds.txt`); regenerates the 83 MB witness file in the temporary directory; B&B for p = 25..36 | 34 s | `ladder.out`, `PASS_EXPQJ_LADDER` |
| `qj/`: `bash ladder.sh --full` | the same plus p = 37..41 | about 12 min more (p = 41 alone ≈ 11 min) | — |
| `qj/`: `uv run -q --no-project --python 3.12 --with numpy python -I -B check.py` | independent numpy checker (shares no code with the C programs): thresholds root 1/25/37, D7 14/42, others ≥ 51; B&B for p = 25..36 with identical Δ and node counts | 9.4 s | `check.out`, `PASS_EXPQJ_CHECK` (also prints `CHECK PASS`) |
| `r3/`: `uv run -q --no-project --python 3.12 --with numpy python -I -B csys_check.py 1` | c-system extra dimension on structured and random supports; controls | 11 s | `csys_check.out`, `PASS_EXPR3_CSYS` |
| `r3/`: `uv run -q --no-project --python 3.12 --with numpy python -I -B members.py 1 9` | row-wise builder, exact κ-members, Lemma 2, label credit ≤ 9, mincost ≥ 72 | 131 s | `members.out`, `PASS_EXPR3_MEMBERS` |
| `nz/`: `uv run -q --no-project --python 3.12 --with ortools python -I -B cluster.py 900` | cl(P120) = 14, CP-SAT optimal (one worker) | 46 s | `cluster.out`, `PASS_EXPNZ_CLUSTER` |
| `nz/`: `python3 -I -B certhc.py 1` | the F2 packing certificate reaches 72 = cost on the P120 frame decoder | 86 s | `certhc.out`, `PASS_EXPNZ_CERTHC` |

Library files: `qj/e8m.py` (self-contained label model, checked against B_rs = (1 + ⟨r,s⟩) mod 2), `qj/mkin.py`, `qj/mkall.py`,
`qj/thr.py`, `qj/pairu.c` (pair witness unions: 25 same line, 26 orthogonal pair, 30 odd pair; usage
`./pairu w67.txt sa jb sb pmax`); `r3/r3lib.py`, `r3/credit.py`, `r3/orb4.json` and the expSM modules (`e8.py`, `kfam.py`,
`maxdef.py`, `cliques_lib.py`, `orbits.py`, `catalog.py`, `subsys.py`); `nz/e8.py`, `*/e8map.json`.

**Large files are not shipped.** The witness file `w67.txt` (83 MB) is regenerated by `ladder.sh`:
```
gcc -O2 -o cand cand.c && python3 -I -B mkall.py allroots_in.txt
./cand allroots_in.txt 41 | awk '$2==6||$2==7' | sort -u > w67.txt    # 1,943,160 lines
```
Not included: the bi-frame annealers (`ann.c`, `annL.c`, `gen*.py`) and their decoder files; expR3's residual searches
(`badsearch*.py`, `rtangent.py`, `resid*.py`); expNZ's ratio anneals. Their results are quoted above.

Sources: expQJ, expR3, expNZ. Prepared with Claude Code (Claude Opus 5.5) assistance. Not independently reviewed.

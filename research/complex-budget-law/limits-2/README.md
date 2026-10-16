# Limits, round 2: Conjecture Θ partial proof, frame routes closed on exact instances, and the remaining lever

Status: LOCAL-ONLY research package, author-checked, not independently reviewed. **It makes no κ claim.** Base: main
`0fdd7d82`. Multipliers are over κ = 0.0017 (our #2679 is at 0.0017149; the selected value is 0.0016727). This round
continues `../limits-1/` (read that index first: the single-letter law, Problem X on each side, Conjecture Θ, the bootstrap
ceiling and the supplier split are defined there). Like limits-1 it is a set of scoped closures and open problems, not a
supplier proposal.

## The picture in one paragraph

Round 1 reduced κ ≥ 0.1 to Problem X on each side and left Conjecture Θ, the bit-side schedule and frames as the open items.
Round 2 attacked each one. **Θ** is now proved to second order on the whole flat manifold (Theorem L, an inertia count) and
exactly on a new global class (Theorem P), and reduced globally to a one-parameter monotonicity statement (Conjecture M) that
survives every adversary; the natural fallback (the linear Schur gauge) is false. **Exact instances** exist: the odd-q
Hermitian code on PG(2,9) (Braverman–He's 166-session code) has an exact helper-free schedule, already published on #2838; at
q = 25 and 49 the best schedules miss by +31 and +157. **Frames** fail on that exact instance under every screen (345 GRS-F₂
ranks, 41 two-class partitions, 184 rational rays, all null), and the instance is an F₉ circuit that the binary contract cannot
take (F₉ has additive exponent 3). **The contract's own premises** were audited: linearity and exactness are forced
(Lemmas N1, B1–B3), dropping the MDS clause closes at ×1.0, mixed local rates are not a lever, and an exact char-0 schedule is its
own realisation, so it already forces s ≥ Wm for rank-metric frames. **The complex side numerically** has a certified Boot
window above E/6 from q = 7 (and above the a = 1/9 threshold from q = 19), but non-equivariant complex MDS codes fall away from
it (r/E 0.154 → 0.091 → 0.066 at q = 3, 5, 7). The one lever left standing is **non-reversible networks (erasure priced as a
local pass), which admit full odd-degree planes**: full PG(2,16) has a bit-side dimension of s = 0.19392 > 0.193 (×60.6
conditional). Its first test case is closed (`full-pg16/`): the best schedule has \|R\| = 909, an exact one needs \|R\| = 900
because dead-session completion fails for this sound code, and the complex Markov cap on that host is 0.18147 < 0.193. The
lever survives only on larger planes or with non-flat local codes.

| Lane | Verdict | Realised | Ceiling left | Kind of evidence |
|---|---|---|---|---|
| [`theta-partial/`](theta-partial/) | **Conjecture Θ: partial proof.** Theorem L (second-order Θ at every flat code, any even degrees, via n₊(K) ≤ E − 1 from Lemmas A, B); global reduction to Conjecture M (ray monotonicity, open, stop test negative); **Theorem P** (exact Θ̂₃ criterion at the point-flat gauge) with **P1/P2** (a global class where Θ holds); the linear-gauge fallback **refuted** (0.1253568 > 1/8 at ‖δ‖_rms = 0.6 on a Γ ≻ 0 base) | ×0 | Θ is a ceiling statement: proving it closes girth-8 complex codes at s ≤ 1/8; refuting it needs a ray with Φ′ > 0 | theorems, exact moments, adversarial numerics |
| [`exact-instances/`](exact-instances/) | **q = 9 EXACT** (published, [#2838](https://github.com/CrocSwap/integer-mult-bounds/discussions/2838)); q = 25 excess **+31**, q = 49 **+157** (relative excess 0.9% → 0.6%); frame screens on the exact q = 9 instance **all null** (345 GRS ranks, 41 two-class partitions, 184 rays); contract mismatch: F₉ has exponent 3, no binary encoding | ×0 | conditional ×37 at q = 9 (needs frames, which the screens exclude); frame route ×1.0 | exact construction, exact ranks and gcds |
| frame-trace (expXS) | **already in [`../limits-1/suppliers/`](../limits-1/suppliers/)**: Lemma T (rank-≤1 frame segments summing to I are trace-1 idempotents), so the exact MC₅ pattern at m = 5 is infeasible on PG(2,16)−M; robust bound s ≥ W − \|X\| + MC_{m−1}. Not repeated here | ×0 | frames with rank-2 segments | lemma, certificate |
| [`contract-slack/`](contract-slack/) | linearity forced (**Lemma N1**), exactness forced (**Lemmas B1–B3**), translation gates frame-compatible (**Lemma G**); erasure is priced as a local pass and mixed rates are not a lever (PG(2,q ≤ 9), expXW); **dropping the MDS clause closes at ×1.0** on PG(2,5), PG(2,7) (every hit above E/6 is a degeneration); an exact schedule with char-0 gates gives s ≥ Wm by its own realisation; **live lead: full PG(2,16) with 8→9/9→8 gates, s = 0.19392 > 0.193** | ×0 | ×60.6 bit-side dimension ceiling at full PG(2,16), conditional on an exact schedule, frames and a complex side; closed for the Hermitian route by `full-pg16/` | lemmas, exact checks, numerics |
| [`complex-search/`](complex-search/) | **spectral Boot bound** Boot(PG(2,q)) ≥ (q+1−√q)(q²+q+1)/4: room above E/6 for q ≥ 7, above the a = 1/9 threshold for q ≥ 19; complex MDS ladders: healthy r ≈ 17 on PG(2,5), ≈ 30 on PG(2,7); r/E falls with q | ×0 | PG(2,q ≥ 19) field-free window (×142 as q → ∞), unreached by any code | theorem, numerics |
| [`full-pg16/`](full-pg16/) | **closed for the Hermitian route.** Best schedule \|R\| = **909** (Singer cyclic m = 11, c = 91; gate-by-gate F₁₆, rank(Φ − I) = 9); the target \|R\| ≤ 906 was wrong: Φ − I has full support (minimum dead cover 906, 3 live sessions), so an exact schedule needs \|R\| = 900; Boot ∈ [900, 901]; new cost identity for unequal k; complex Markov cap at traces 8/17, 9/17 **M = 0.18147 < 0.193** (κ ≤ 0.060 on this host) | ×0 | ×35 at most (exact F = 900 *and* a complex code at the Markov cap) | exact simulation and cover, exact merge DP, heuristic search, LP |
| [`pending/`](pending/) | placeholder: the record column-operation ILP (still running) | — | — | — |

**Realised multiplier: ×0 for every lane.** Nothing here moves κ.

## What is proved

1. **Θ at second order everywhere on the flat manifold** (`theta-partial/` Theorem L, Lemmas A, B): at every two-sided-flat MDS
   code (girth ≥ 8, any even degrees) the second-order code/gauge max–min is ≤ 0; the bound n₊(H₀) ≤ E is sharp
   (cos 2θ = −1/9). The refutation route of round 1 ("flat code with Γ ⪰ 0 and a positive normal Schur direction") is closed
   wherever H₀ is nonsingular.
2. **Global reduction** of Θ to flat codes plus a relative gauge δ, and **Theorem P**: an exact criterion for Θ̂₃ ≤ 1/8 at the
   point-flat gauge, with **P1** (no point on two non-flat lines ⇒ Θ̂₃ = 1/8) and **P2** (two per point with opposite
   deviations) as a global class where Θ holds, strictly larger than round 1's two-sided-flat class.
3. **The q = 9 odd-characteristic Hermitian code is exactly schedulable** (`exact-instances/`; gate-by-gate F₉ simulation,
   rank(Φ − I) = 0), and **S(25) = 3449** is confirmed independently by simulating the cyclic schedule (rank(Φ − I) = 31).
4. **Contract facts**: no binary encoding of F₉-linear gates (exponent 3 vs 2); a sound exact char-0 schedule is its own
   monomial realisation (s ≥ Wm for every rank-metric frame); at girth 6 the flow screen F_* < W is automatic once S > E/6.
5. **Forced premises** (`contract-slack/`): gates must be linear under any non-monomial frame (N1); border and approximate
   frames give nothing in the additive contract and cost ≥ ×2 precision per level in the complex supplier (B1–B3); the
   commutant of the phase frames is the translation algebra (G); zero scratch and erasure cost 0 in the exponent, and the
   Problem-X threshold does not depend on local rates (R1–R5).
6. **Spectral Boot bound** for every PG(2,q) (`complex-search/`), field-free.
7. **Full PG(2,16) Hermitian code** (`full-pg16/`): the Singer schedule with |R| = 909 realises dim W = 900 (rank(Φ − I) = 9,
   gate-by-gate over F₁₆), Φ − I has full support and the minimum dead cover is 906, so dead-session completion fails and
   an exact schedule needs |R| = 900; the cost identity cost(σ) = ½(ascent + #rising lines − #rising points) for unequal k;
   Boot ∈ [900, 901]; Markov cap 0.18147 at traces 8/17, 9/17 (under expRF's flat-leverage premise).

## What is open

| Item | Side | Exact question | Evidence so far | Stop test |
|---|---|---|---|---|
| **Conjecture Θ** | complex | every MDS code at girth ≡ 0 mod 4 has a gauge with Θ̂_{g/2−1} ≤ 1/g | proved at second order (Thm L) and on the P1/P2 class (Thm P); one-sided gauges fail, so interior gauges are needed | a code with min_s Θ̂₃(x_s) > 1/8 on the point-flat → line-flat segment while Θ̂* ≤ 1/8 |
| **Conjecture M** | complex | Φ(t) = inf_x Θ̂₃(C₀; e^x, e^{x+tδ}) is non-increasing (M ⇒ Θ) | 39 rays monotone; 9 adversary runs on J = Φ(δ) − Φ(0.8δ) all negative; no per-vertex identity exists | a precisely re-solved ray with Φ′ > 0 |
| **Problem X(ℂ, 1/9)** | complex | a schedulable MDS complex Tanner code with s ≥ 0.193 (g = 6) | field-free window exists on PG(2,q ≥ 19); codes reach only r/E ≈ 0.066 at q = 7 and fall with q | reopen only for a structured non-flat family with r ≥ 0.12·E on PG(2,7) |
| **Exact schedules at q = 25, 49** | arithmetic (MSP) | is Boot = S, with a feasible transportation? | +31 / +157 schedulable, Boot ≤ S + 16 / S + 81; cyclic is a strict local optimum of F | q = 25: F-aware search ≥ 1 h staying ≥ 3480 |
| **Frames on exact instances** | arithmetic (MSP) | a frame family with positive nullity on the exact q = 9 instance | all screens null; GL(m, F₂) address-map frames and universal certificates not covered; char-2 realisability of the incidence undecided | CryptoMiniSat on the exact incidence, 1 h |
| Non-reversible odd-degree planes beyond PG(2,16) | bit and complex | an exact schedule (F = dim) for a full odd-degree plane with dim/E > 0.193, *and* a complex code above 0.193 on the same host | PG(2,16) closed for the Hermitian route (`full-pg16/`: best \|R\| = 909 ≠ 900; Markov cap 0.18147 for flat-leverage codes); non-flat local codes and larger planes untested | first a paper bound for non-flat local codes at unequal traces; schedule search only if it exceeds 0.193 |
| Problem X(F_{2^e}, 1/9) | bit | schedule for PG(2,16)−M / PG(2,64)−M (round 1) | unchanged this round | `../limits-1/bootstrap/` |
| Translation-term gates | both | does one X_w term per gate lower the exact minimum s on the control instances? | not run (P ≈ 2%) | minima equal on all control instances |

## The remaining lever

Of the contract's unforced premises (expL8), only **reversibility** survives: zero-born scratch and free erasure cost nothing in
the exponent, keep every screen, and let non-square gates use full odd-degree planes. That is what moves full PG(2,16) above
the a = 1/9 line on the bit side (dimension 900/4641 = 0.19392, conditional κ ≤ 0.103, ×60.6). On that host it fails twice
(`full-pg16/`): the schedule must be exact, |R| = 900 rather than ≤ 906, because dead-session completion cannot localise the
defect of a sound code (best 909, the Singer family is exactly optimal at 909); and the complex side cannot reach 0.193 there
with flat-leverage codes (Markov cap 0.18147, a ≤ 0.064, ×35.4, even with a perfect schedule). So the lever is left only for
larger planes (PG(2,q ≥ 19), where the field-free complex window exists) or for non-flat local codes, and every route in this
package ends at the same open item: a complex code above 1/g. Conjecture Θ would rule that out at girth 8; at girth 6 (the
planes) no proved bound rules it out, and no code found comes close.

## Replay

```
bash run_all.sh            # everything (19 checks), about 5 min on one worker
bash run_all.sh --quick    # skips the 4 runs longer than about 20 s (15 checks, under 1 min)
```

Stdlib scripts run on system `python3`. Scripts that need numpy, scipy, networkx or ortools run through
`uv run -q --no-project --python 3.12 --with …`; each lane README gives the exact command. All scripts locate their inputs
relative to their own path (scripts that import a sibling insert their own directory, since `python -I` omits it). Inside this
branch, `theta-partial/xr` imports expXK/expXG from `../limits-1/complex-dual/{xk,xg}/`, and `contract-slack/xv`,
`complex-search/xp` read `../limits-1/complex-dual/xg/graphs_small.json`. All checks were run on one worker under a 4G memory
cap (`systemd-run` scope); BLAS/OpenMP threads are pinned to 1.

## Sources

Condensed from local scouting runs: expXR (theta-partial); expXQ, expXU (exact-instances); expXS (frame-trace, in limits-1);
expL8, expXW, expXV (contract-slack); expXP (complex-search); expXX (full-pg16); the column-operation ILP (pending, still running).
The q = 9 schedule artifact is the one published on #2838. Upstream files are not used, except that `exact-instances/xu/contract.py`
(not replayed) calls upstream `scripts/finite_bit_contract.py`.

Prepared with Claude Code (Claude Opus 5.5) assistance, using many sub-agents. Not independently reviewed.

## Run record

Full replay on 2026-10-16: `CAP=4G capped bash run_all.sh` (one `systemd-run` scope, 4G memory cap, one job at a time, one
worker): **19/19 checks PASS, exit code 0, 291 s wall time.** Quick replay (`--quick`, skipping `check_xr.py`, `check_xr2.py`,
`cyc.py 49` and the PG(2,5) `cmds.py` rung): 15/15 PASS, 27 s. Individual runtimes vary with load; the lane README tables give
values from isolated runs. Replayed outputs that have a scouting counterpart were compared with it: `check_xr`, `check_xr2`,
`check_l8`, `markov17` and expXQ's q = 9 / q = 25 verifier lines are identical; `schedule9.json`, `cyc{9,25,49}.txt`,
`graph{9,25,49}.txt` regenerate byte for byte; `sim16.py` reproduces `sched909.json` (order, arcs, R, rank, cover size).
Two differences, both without effect on a verdict: `xw.py`'s random mixed-profile draws were salted per process in scouting
and are now seeded deterministically (column C only); the PG(2,5) `cmds.py` rung replays to f = 5.0e-7 against the logged
1.0e-7. Not part of the replay: the scouting searches each lane README lists as "included but not replayed".

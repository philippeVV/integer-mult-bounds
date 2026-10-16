# Scoped closures, round 3 (complex budget law)

Status: LOCAL-ONLY research package, author-checked, not independently reviewed. **It makes no κ claim.** Base: main
`0fdd7d82`. Multipliers are over κ = 0.0017. Guidance refreshed against main `fbdbcb53` (priority revision 5,
`docs/research/research-family-status.md`). The claims on the board are at 0.0017187 (#2833/#2843), and the selected value is
0.0016727.

Each lane closes, caps or bounds one structural route to a κ multiplier, except `char2-tanner-family/`, which is a rev-5
**family record** (an expansion decision), not a closure. Conventions are those of `../closures-2/`. Every lane directory holds:
- a README with the verdict, the changed hypothesis, a proof sketch, the multiplier and the stop test;
- minimal scripts, which locate their inputs relative to their own path (repository files through `IMB_REPO` or `parents[4]`);
- their outputs (`*.out`). Every script prints a `PASS_` marker.

| Lane | Verdict | Realised | Ceiling left | Kind of evidence |
|---|---|---|---|---|
| [`peel-dead/`](peel-dead/) | **closed**: the u = 1 peel is dead (every single climb is exactly +I, Theorem B), and so is the peel with any **fixed** number of extra climbs (δ ≤ C(k)R, Theorem F). **≤ W + 2 climbs ⇒ s ≥ Wm** (Theorem G); the δ = R hyperbola is not realisable; with **any number of independent** extra climbs, locality and single use give δ ≤ 0, so a = 0 (Theorems L, U). Small two-port gadgets: routability conjecture open, scoped, no κ impact | ×0 | dependent extras (k ≥ 3) and strict deficits outside the level setting are open | theorems, exhaustive small searches, exact MITM |
| [`large-ports-trace-bound/`](large-ports-trace-bound/) | **closed (idea 1 at ×1.0)**: rank N ≥ uV/h over any char-2 field (Theorem A); PSD local ϑ bound (Theorem B); E8-40 (factor 5/4) has no useful realisation; Q/F2 separation σ ≥ 1 numerically on every field-line target. **Field-line / polar large ports retired** | ×1.00 | reopen only with non-field labels (PR ∈ 1..u − 1, density ≥ 0.1); the real non-symmetric case is open exactly at 2-adic free cycles | theorems, exact certificates, SDP, numerics |
| [`laminar-decoders/`](laminar-decoders/) | **closed**: no laminar (helper-free, direct-gate) decoder compresses. Bowtie obstruction; rank lemma; span MILP proves the identity optimal in the compressing span | ×1.00 | ≤ ×1.05 side effect (partial laminarity on today's word, uncharged) | exhaustive certificate, MILP, pricing |
| [`endpoint-labels-both-banks/`](endpoint-labels-both-banks/) | **closed** in template form (rank invariant; oracle strictly worse with any Y label) and in one-block form at every dense decoder (Theorem W3: threshold back to γ < ½) | ×1.0 | one unconstructed object: two consecutive decoders with ≥ 1/7 flat rows at cost < (2/3)V | theorems, exact oracle, exhaustive toy search |
| [`char2-tanner-family/`](char2-tanner-family/) | **family record, decision PAUSE.** BIT-ONLY, so ×1.0 in κ on its own. Every tested frame family certified ×0 (GRS, ray test, F1 dead); no exact schedule (784 vs 772; divisibility, cut + viol); no code beyond 772. expQR and expRA filled in (MC_m ≥ Wm for m ≤ 4 not certified; MC_5 leaves room) | ×0 | bit-side ×10–25 conditional (worth ×1.0 in κ until a complex counterpart exists) | exact certificates, simulations, mod-p ranks |
| [`complex-coding-seeds/`](complex-coding-seeds/) | **closed (expL6 C1)**: no char-0-only Tanner gap above threshold for flat-leverage (ζ / Walsh–Hadamard) codes. Chebotarev rigidity in char 0 (the char-2 gap is a mod-2 Chebotarev failure); T1: S0/E ≤ 1/4 on any layout; T3: a ≤ 0 at girth 4, 8, 12, 16; every char-0 gap reappears mod 2. **The complex side has no coding-seed route** | ×0 | girth-6 sliver (S/E ≤ 2/11, a ≤ 0.065) and girth 10 (a ≤ 0.017), both unreached | theorems, exact Fourier/mod-p ranks, LP |
| [`small-units-lp/`](small-units-lp/) | **update to closures-2/small-units-h5, now PAUSED**: trace rows are invalid (spreading counterexample); exact idempotent caps; tight-regime floors 5/16 (h = 6) and 1165/2837 (h = 5) fail the stop lines | ×1.00 | unchanged: ×3.03 ('best') at h_c = 5; ×16.7 / ×14.7 planning ceilings for rank-additive units | exact linear algebra, rational LP duals |
| [`decoder-law-subsets/`](decoder-law-subsets/) | **partial closure, cross-posted to the decoder thread #1587 (extends PR #968)**: 8-term decoders have γ ≥ ½ on every \|Q\| ≥ 79 (p ≤ 41; independent checker p ≤ 36); Theorem R3 closes rank-9 lifts on supports ≤ 26 lines; L1–L3 tools. P104 (56/104) stays the cheapest ratio | ×0 | p ≥ 42 open (span-only relaxations provably fail for \|Q\| ≤ 47); R3″ residue; any decoder at \|Q\| ≤ 95 needs gain decay | exact census in C plus independent numpy checker, proofs |

## Replay

```
bash run_all.sh            # everything (67 checks), about 11 min on one worker
bash run_all.sh --quick    # skips the 5 runs longer than about 1 min (about 5 min)
```

Stdlib scripts run on system `python3` (3.14 here). Scripts that need numpy, scipy, ortools, networkx, mpmath, cvxpy or
python-flint run through `uv run -q --no-project --python 3.12 --with …`; each lane README gives the exact command. C programs
(`endpoint-labels-both-banks/bfs.c`, `decoder-law-subsets/qj/{sub,cand,bb2}.c`) are compiled with `gcc -O2` by the lane's shell
driver into a temporary directory. The 83 MB witness file of `decoder-law-subsets/qj/` is regenerated, not shipped.

All checks were run on one worker under a 4G memory cap (`systemd-run` scope); the runtimes are in each lane's README. Results
of the full run: see "Run record" below.

## Late inputs (folded in)

- **expQX**, **expQY**, **expRG** and **expRH** have reported and are folded into `peel-dead/`; **expQA** into `large-ports-trace-bound/`;
  **expQR** into `char2-tanner-family/`.
- **expRA (MC_4)** has reported and is filled into `char2-tanner-family/` (not certified; the decision stays PAUSE).

## Sources

The lanes were condensed from local scouting runs: expQQ, expQH, expQX, expQY, expRG, expRH (peel); expQL, expQA, with expQF and expNM as
cross-references (large ports); expQS (laminar); expQT, expQW (both banks); expNC, expNI, expNJ, expNK, expNY, expQB, expQD,
expQE, expQI, expQK, expQR, expQU, expQV, expRA, expL6 (Tanner family); expRB, expRF (complex coding seeds); expNG (small units); expQJ, expR3, expNZ (decoder law).
Inputs are copied into each lane where needed (pickles converted to JSON). Upstream files are used in place and unchanged:
`research/biframe-decoder-ladder/RESULTS.json`. Within this branch, lanes import `../../framework-ceiling/` (through the
closures-2 copy of `kappa.py`), and `char2-tanner-family/l6/estimates.py` imports `../../peel-dead/hyp_drop.py`.

Prepared with Claude Code (Claude Opus 5.5) assistance, using many sub-agents. Not independently reviewed.

## Run record

Full replay on 2026-10-15: `CAP=4G capped bash run_all.sh` (one `systemd-run` scope, 4G memory cap, one job at a time, CP-SAT at one worker),
**67/67 checks PASS, exit code 0, 634 s wall time.** Individual runtimes vary by up to about 2× with machine load; the per-lane
README tables give the values from the isolated runs. Not part of the replay: `decoder-law-subsets/qj/ladder.sh --full`
(p = 37..41, about 12 min; p = 41 was run once separately, max Δ = 6 ≤ 8), and the scouting searches each README lists as not
included.

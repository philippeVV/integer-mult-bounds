# Where the limit is: single-letter theory, suppliers, and the remaining open problem

Status: LOCAL-ONLY research package, author-checked, not independently reviewed. **It makes no κ claim.** Base: main
`0fdd7d82`. Multipliers are over κ = 0.0017 (our #2679 is at 0.0017149; the selected value is 0.0016727). Guidance refreshed
against main `fbdbcb53` (priority revision 5, `docs/research/research-family-status.md`): this package is a set of scoped
closures and one exact open problem, recorded as a family decision ("new coding, geometry and recursive boundaries"), not a
supplier proposal.

## The picture in one paragraph

The recursion is not the bottleneck. κ is decided by **the best finite seed on each side**: the varying-seed law makes κ a
log-Cesàro mean of seed exponents, tensor products and mixtures cannot beat their best member, constant overheads amortise,
and the repo recursion already amplifies a coding gap 10–20× better than the published graph-product amplifiers. One seed with
a ≥ 1/9 on each side would give κ ≈ 0.1 (×59). Integer multiplication cannot drop either side: the movement ("bit") supplier is
char 2 in every recipe, and the arithmetic ("complex") supplier cannot have exponent 2 (Lemma A); every outer whose
characteristic grows is char-0 class. So the problem is **Problem X on each side** — a schedulable MDS Tanner code with excess
s = dim C/|E| above 1/(1+(g−1)^{1−a}) — and schedulability itself forces dim C = Boot(Γ), a field-free combinatorial ceiling.
On the bit side the dimension exists (Hermitian PG(2,64)−M, s = 0.2154) and the schedule is open. On the complex side every
proved bound points at a ≤ 0, but the general statement is open: Conjecture Θ at girth 8 (tight at ETF codes) and girths 6/10
with no dual-only route.

| Lane | Verdict | Realised | Ceiling left | Kind of evidence |
|---|---|---|---|---|
| [`single-letter/`](single-letter/) | **theory**: κ = log-Cesàro mean of seed exponents (Theorem 1); one seed with a ≥ 1/9 per side suffices for κ ≥ 0.1; the repo recursion beats BGS by ×9.9–20.7 on the BH seed; **Problem X** stated exactly | ×0 | conditional ×142 (Hermitian, either side alone does nothing), outer-capped ×147 | theorem, simulation |
| [`complex-dual/`](complex-dual/) | **ceiling map**: gauge-invariant dual Θ̂*; girth 4 closed for all codes; flat codes closed at every even girth ≤ 16 (s ≤ 1/g); two-sided-flat codes closed at girth 6 and 8; **Conjecture Θ (girth 8) open and tight**; girths 6, 10 open for non-flat codes, where the moment dual is provably insufficient | ×0 | the open cells of the map (non-flat codes at g = 6, 8, 10, …) | theorems, exact identities, adversarial numerics |
| [`bootstrap/`](bootstrap/) | **theorem, any field**: dim C ≤ Boot(Γ); schedulable ⇒ dim C = Boot; X(·, 1/9) excluded on PG(2,q ≤ 13), PG(2,16)−M, W(q ≤ 7) and all degree-4 graphs tested. expXL: Boot(PG(2,16)−M) ∈ [772, 784], undecided; cutwidth cannot decide it; exact cut-ascent identity | ×0 | PG(2,q ≥ 17), W(q ≥ 11): Boot/E up to ~¼ (×142) | theorem, greedy certificates |
| [`suppliers/`](suppliers/) | **closed for "bit side only"** (Lemma A); B′, Greenberg, Ax put growing-characteristic outers in the char-0 class; MSP is the one escape; phase frames reduce to Sym_m(F₂) and Theorem Ψ forbids savings at m ≤ 3 in every field; odd-q Hermitian S(q) closed form (166, 3449, 26434). expXN: odd-p phase frames give nullity 0 at m = 4–8, closed at ×1.0; expXS: rank-≤1 frames at m = 5 infeasible on PG(2,16)−M (trace lemma) | ×0 | ≤ ×1.01 with today's bit side; ×97–187 only if both sides get frames and schedules | lemmas, exhaustive checks, exact ranks |
| [`lattices/`](lattices/) | **closed (decoder lane)**: Barnes–Wall invalid by the valuation lemma; Lorentzian E_n gives γ_n = 6/(n−1) → 0 but γ·(h−1)/u ≡ 6 > E8's 4.57; fixed windows never exceed E8's 120 lines; Hermitian ports are a tensor copy of real ones | ×1.0 | ×2.80 selected / ×5.02 data floor (u = 1 free-decoder row) | exact checks, CP-SAT, pricing |
| [`seeds-negative/`](seeds-negative/) | **closed**: hexagons H(3)/H(4)−M/H(5) (S/E ≤ 0.097 vs 0.232 needed); Singer-equivariant complex codes on PG(2,17) (≤ 81 exact, 120 numerical vs 1005; equivariant-degeneracy lemma); expander ceilings (a ≥ ½ ⇔ S/E ≥ 1/(1+√(g−1)); cut/flow cannot close it, codes decide). expRB/expRF: see `../closures-3/complex-coding-seeds/` | ×0 | conditional only (×14–207), all needing codes that are not known | exact ranks, numerics, LP |

**Realised multiplier: ×0 for every lane** (×1.0 for the lattice lane, which only changes decoders). Nothing here moves κ.

## What is proved

1. **Single-letterisation** (`single-letter/`, `complex-dual/` Thm E): Ψ and Θ̂ are multiplicative under tensor products and the
   varying-seed recursion gives the log-Cesàro mean, so a seedwise bound closes all tensor products and mixtures.
2. **Two suppliers are unavoidable for integers** (`suppliers/` Lemma A, B′; B, C modulo Greenberg and Ax): κ ≤ min(a_mov,
   a_arith), and a_arith is char-0 class for every outer except multi-small-prime.
3. **Bootstrap ceiling in every field, and schedulable ⇒ dim C = Boot** (`bootstrap/`).
4. **Complex ceilings** (`complex-dual/`): gauge invariance is necessary (Lemma D); s ≤ ¼ for all codes (girth 4 closed);
   s ≤ 1/g for flat codes at g = 4…16; two-sided-flat codes give Θ̂ = 1/6 at girth 6 and 1/8 at girth 8 (Thm F, Thm 2, with the
   exact m₃ identity).
5. **Phase frames** (`suppliers/`): they are polar forms in Sym_m(F₂); the cost is char-blind (Lemma Y); GRS-F₂ forces s ≥ Wm
   on char-2-realisable incidences; Theorem Ψ forces s ≥ Wm at m ≤ 3 on every realisable circuit.
6. **Decoders** (`lattices/`): valuation lemma; tensor law for one-valuation Hermitian ports.

## What is open (the remaining problem, by side)

| Item | Side | Exact question | Evidence so far | Stop test |
|---|---|---|---|---|
| **Problem X(F_{2^e}, 1/9)** | bit | does a Hermitian Singer code (PG(2,16)−M at 772, or PG(2,64)−M at s = 0.2154) have a flow-feasible bootstrap-optimal order, i.e. Boot = dim C? | dimension exists; Boot(PG(2,16)−M) ≤ 784; greedy orders on PG(2,3), PG(2,9) are not flow-feasible; expXL: Boot ∈ [772, 784], PG(2,64)−M Boot/E ∈ [0.2154, 0.2163] | a certified lower bound Boot > 772 closes q = 16 |
| **Problem X(ℂ, 1/9)** | complex | a schedulable MDS complex Tanner code with s ≥ 1/(1+(g−1)^{8/9}) (0.193 at g = 6, 0.151 at 8, 0.124 at 10) | every proved bound gives a ≤ 0; nothing constructed exceeds 1/g (best 10/64 at girth 6, ETF codes at exactly 1/8 at girth 8); only graphs with Boot/E > threshold can host it (PG(2,q ≥ 17), W(q ≥ 11)) | see Conjecture Θ; any instance with s > 1/g is decisive |
| **Conjecture Θ** | complex | g ≡ 0 mod 4: every MDS code has a gauge with Θ̂_{g/2−1} ≤ 1/g | tight at SIC / icosahedral ETF codes; second-order local maximality; W(5) adversary 0.12499917; the complementary gauge alone is provably insufficient | a flat degree-4 code with Γ ⪰ 0 and a positive normal Schur eigenvalue refutes it |
| **Frames for odd-p codes** | arithmetic (MSP) | a phase-frame family (m ≥ 4, non-commuting β) with GRS-F₂ nullity ≥ 0.02·Wm on a non-char-2-realisable incidence | Ψ forbids m ≤ 3; expXN screening m = 4–8 | nullity 0 for all families at m ≤ 8 closes the finite-field route at ×1.0 |

A κ ≥ 0.1 construction needs **both** a bit-side and a complex-side (or arithmetic-side) solution; by `single-letter/` it would
also be a new network-coding theorem (undirected coding gaps ≥ (lg N)^{0.1}).

## Rev-5 family decision

- **Decision: PAUSE all seed searches; RETAIN Problem X as the single open target.** The lanes here close every concrete seed
  family tried (frequency, WH, hexagon, Singer-equivariant, lattice decoders, growing families) and every recipe-level escape
  except MSP. Another seed search is not a different mechanism (rev-5: "not another seed").
- **Reopen evidence:** (i) a proof or refutation of Conjecture Θ; (ii) a complex code with s > 1/g on any graph with
  Boot/E > 1/g; (iii) expXL showing Boot(PG(2,16)−M) = 772 with a feasible order; (iv) expXN finding a frame family with
  positive GRS-F₂ nullity.
- **Obligation ledger** for the one surviving supplier route (MSP) is in `suppliers/`; unknowns there are not zero.

## Replay

```
bash run_all.sh            # everything (38 checks), about 3.5 min on one worker
bash run_all.sh --quick    # skips the 2 runs longer than about 20 s (36 checks, about 2 min)
```

Stdlib scripts run on system `python3` (3.14 here). Scripts that need numpy, scipy, numba, networkx, ortools, cvxpy or clarabel
run through `uv run -q --no-project --python 3.12 --with …`; each lane README gives the exact command. All scripts locate their
inputs relative to their own path (every script that imports a sibling inserts its own directory, since `python -I` omits it).
Inside this branch, `lattices/xa/price.py` imports `../../framework-ceiling/` and two closures-3 copies (`kappa.py`,
`hyp_drop.py`). All checks were run on one worker under a 4G memory cap (`systemd-run` scope); BLAS/OpenMP/numba threads are
pinned to 1 and CP-SAT uses one worker.

## Sources

Condensed from local scouting runs: expXB (single-letter); expXG, expXK (complex dual); expXE, expXL (bootstrap);
expXF, expXJ, expXN, expXS (suppliers); expXA, expXH (lattices); expXD, expXM, expRJ, with expRB/expRF cross-referenced to
closures-3 (seeds). Inputs copied into the lanes: expXE's PG(2,3) r = 8 solution and graph files, expN1's BH incidence, expNC's
`gf.py`. Upstream files are not used.

Prepared with Claude Code (Claude Opus 5.5) assistance, using many sub-agents. Not independently reviewed.

## Run record

Full replay on 2026-10-16: `CAP=4G capped bash run_all.sh` (one `systemd-run` scope, 4G memory cap, one job at a time, one
worker): **38/38 checks PASS, exit code 0, 199 s wall time.** Quick replay (`--quick`, skipping `odd_singer.py 7` and
`h5.py 5`): 36/36 PASS, 123 s. Individual runtimes vary with machine load; the per-lane README tables give the values from the
isolated runs. Every replayed output that has a scouting counterpart was compared with it and is identical (`xb`, `check_xg`,
`check_xk`, `xf`, `psi_check`, the expXA outputs, `part2`, `price2`, the expRJ summary), apart from the appended `PASS_` lines.
One replay finding: expXE's W(7) greedy value (483, Boot/E = 0.1509) was above the exact κ ≥ 0.1 threshold 0.15063; with more
restarts the greedy reaches 479 (0.1497), so the exclusion stands (`bootstrap/README.md`). Not part of the replay: the
scouting searches each lane README lists as "included but not replayed".

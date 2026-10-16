# Scoped closures, round 2 (complex budget law)

Status: LOCAL-ONLY research package, author-checked, not independently reviewed. **It makes no κ claim.** Base: main
`0fdd7d82`. Multipliers are over κ = 0.0017. Guidance refreshed against main `fbdbcb53` (priority revision 5). The claims on the board are at 0.0017178 (#2802/#2818), and the selected value is
0.0016727.

Each lane closes, caps or bounds one structural route to a κ multiplier. Every lane directory holds:
- a README with the verdict, the changed hypothesis, a proof sketch, the multiplier and the stop test;
- minimal scripts, which locate the repository relative to their own path (override with `IMB_REPO`);
- their outputs (`*.out`).

Every script prints a `PASS_` marker.

| Lane | Verdict | Realised | Ceiling left | Kind of evidence |
|---|---|---|---|---|
| [`flat-tower/`](flat-tower/) | **closed**: n ≤ 4 pairs per replica. A shared stage needs a rank-8 repair; 2n+1-call words have only rank-1 blocks. **Corrects the earlier ×1.32 row.** | ×1.00 | ≤ ×1.05 (an unknown parallel word with e ≥ 2) | exact GF(2) on the real routes, a theorem, pinned engines |
| [`zll-recursion/`](zll-recursion/) | **closed**: recursion, tensor powers and cross-level overhead sharing are ×1.000. Ψ_{S⊗S′} = Ψ_S Ψ_{S′}; Lemma W′ and the switch-cycle bound; the Kronecker overlap lemma; no carrier/source overlap on #1094 or #589. | ×1.00 | generous residual ×1.016; ideal ×2.77 (centre) and ×9.65 (centre + helpers), both unreachable | proofs plus exact checks against the upstream contract scorer |
| [`sibling-helpers/`](sibling-helpers/) | **capped**: joint use of helper registers by sibling calls needs a rate-½ mixed code. Lemmas P and S; the reader cuts are tight but not violated. | ×1.00 | ×1.41 (k = 2, perfect code); ×1.33 with junk charged | pinned engines, exact census mod 2⁶¹−1 |
| [`p120-rank8-rigidity/`](p120-rank8-rigidity/) | **closed for line ports at h ≤ 12**: idea 4 (sparse transfers) and the corrector-rank lever. Every rank-≤8 P120 decoder is Λ(G/2)Λ⁻¹ (T = F). Every structured sparse support with d̄ ≤ 30 needs rank ≥ 12. Half-integral decoders have T = F. | ×1.00 | none known | theorem with a DRUP-certified finite step; exact α bounds |
| [`tm-pass-programs/`](tm-pass-programs/) | **closed**: the whole-volume masked-lag pass model is κ-blind (needs Θ(u/log k) passes). XOR = route at every small size solved. | ×1.00 | none | proof plus exact SAT tables |
| [`field-line-large-ports/`](field-line-large-ports/) | **stop test fired**: the best exact decoder has γ = 1 in all 8 polar families. The minrank law V/ω is **numerical only**. | ×1.00 | ×24–36 at the unexcluded α floors | exact decoders; LM evidence |
| [`omega-h-floor-false/`](omega-h-floor-false/) | **refutation**: C ≥ ω·h is false (P112: C ≤ 64 < 72; F_7: 50 < 56; and others). The V/ω screen is invalid; use per-object floors plus monotone gain. | ×1.00 | — | exact certificates |
| [`small-units-h5/`](small-units-h5/) | **OPEN**: no unit with x < 1 found, and no certificate of non-existence. | ×1.00 | ×3.03 ('best') at the LP floor; break-even needs x ≤ 0.7275 | pair scoring; exact verifier |

## Replay

```
bash run_all.sh            # everything, about 15 min on one worker
bash run_all.sh --quick    # skips the 8 runs longer than about 1 min (about 4 min)
```

Stdlib scripts run on system `python3` (3.14 here). Scripts that need numpy, ortools, python-sat or python-flint run through
`uv run -q --no-project --python 3.12 --with …`.

All checks passed on the full run (one worker, 4G memory cap). The runtimes are in each lane's README.

## Sources

The lanes were condensed from local scouting runs (expNV, expNL, expQC, expQG, expQM, expNU, expNO, expQF, expNX, expNM,
expQP, expNF). Their inputs are copied into each lane where needed: gzipped networks, the CNF and DRUP proof, and the
candidate schedules converted from pickle to JSON. Upstream research packages are imported in place and unchanged:
- `complex-local-word-ceiling/ceiling.py`;
- `bit-revision3-levers-seed03/levers.py`;
- `scripts/finite_bit_contract.py`;
- the #1094 selected replay;
- the #589 tuned word;
- `direct-delivery504-global` DATA-PHASES.

Prepared with Claude Code (Claude Opus 5.5) assistance, using many sub-agents. Not independently reviewed.

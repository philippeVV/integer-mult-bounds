# Large ports (idea 1): the V/ω law is a theorem in characteristic 2 and for PSD decoders; field-line / polar large ports retired (×1.0)

Status: scoped closure, author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`.
Multipliers are over κ = 0.0017.

## Verdict

| Item | Value |
|---|---|
| Realised | **×1.00**. The stop test fires: no label system at h ≤ 17 combines a realisable super-frame factor > 1.2 with label fraction k ≥ 1/4 |
| Best graph-level factor at ω = 4 | **5/4**, from E8-40 (40 E8 root lines, no 5 pairwise non-orthogonal), certified exactly. Its only working label realisation is the P120 restriction at h = 9, u = 1, where φ_eff = uV/(hK) = 0.56: no gain |
| E8-40, hypothetically realised in PG(3,8) / PG(3,16) | ×2.11 / ×2.06, **below** today's E8 P112 control (×4.05–4.64 in the same formula). Break-even with the control needs φ ≈ 1.45 |
| Theory, characteristic 2 | **closed**: rank N ≥ uV/h for every label system and every decoder over a field of characteristic 2 (Theorem A) |
| Theory, PSD decoders | **closed** for planes (every q), PG(3,4) and P120 (tight at 8); PG(3,8) factor ≤ 1.207 (Theorem B) |
| Open | general non-symmetric real decoders at n ≥ 3. Theorem A pins the gap exactly to decoders with a 2-adic denominator on some free cycle |
| Q/F2 separation screen (expQA) | σ = mr_Q/(θV) ≥ 1 **numerically** on every field-line target at θ ≥ 1/4 (q ≤ 32); rigorous floors **0.33–0.64** at the large targets (PG(2,4) control proved exactly: mr = 7, σ = 21/16). No masked rational matrix below the frame-packing value exists in any screened family |
| **Verdict** | **field-line / polar large ports are retired.** The only reopen route is **non-field labels with PR ∈ {1, …, u − 1} and free-pair density ≥ 0.1** |

This turns closures-2's "the minrank law V/ω is numerical only" (`../../closures-2/field-line-large-ports/`) into a theorem in
characteristic 2 and for PSD decoders, and it records the search that closes idea 1 as a lever.

## Changed hypothesis

expNM's explanation was that field lines have factor 1 because their free density is about 2^{−u}. Density is **not** the
controlling variable (Kneser coordinate labels have density 0.38 and factor ≤ 1 exactly; polarity planes are closed by local
structure). The two invariants that do control it are:
1. the Q-versus-F2 min-rank separation of the free graph (Theorem A);
2. the local ϑ-excess of free neighbourhoods, for PSD decoders (Theorem B).

## Theorem A (char-2 trace bound; exact; any label system)

Labels L_a are u-dimensional subspaces of F2^h with pairing PR(a, b) = rank(L_aᵀ Ω L_b). A decoder has N_aa = 1 and N_ab = 0
on non-free pairs (PR ≠ 0).

**Statement.** Over any field of characteristic 2, rank N ≥ Σ_a PR(a, a)/h = uV/h.

**Proof.** G = (L_aᵀ Ω L_b)_{a,b} has rank ≤ h, so (N ⊗ J_u) ∘ G has rank ≤ h·rank N. Its off-diagonal blocks vanish (N_ab = 0
on non-free pairs, G_ab = 0 on free pairs), so it equals blockdiag(G_aa), of rank Σ PR(a, a). ∎

**Corollaries.**
- (A1) Field lines, char 2: h = uω, so every F2/F_q decoder has rank ≥ V/ω. With frames this is the exact law.
- (A2) Over Q or a number field: if N is gauge-equivalent (N ↦ DND⁻¹) to a 𝔭-integral matrix for some 𝔭 | 2 — exactly when no
  directed cycle of nonzero entries has a product of negative 𝔭-valuation — then reduction mod 𝔭 keeps the mask and the unit
  diagonal, and the bound holds over Q too. Beating uV/h therefore needs an algebraic decoder with a 2-adically negative free
  cycle at every prime above 2. E8's G/2 (entries ±½) is exactly such a decoder.
- (A3) φ_eff := uV/(hK) ≤ mr_F2(F)/mr_R(F): the gain over frames is at most the Q/F2 min-rank separation of the free graph.

## Theorem B (PSD decoders: local ϑ bound)

**Statement.** If N is PSD (a Gram of unit vectors) and θ_a = ϑ(F̄[N(a)]) (the non-free graph induced on a's free
neighbourhood), then rank N ≥ V²/Σ_a(1 + θ_a), and the PSD super-frame factor is ≤ (1 + θ̄)/ω.

**Proof.** rank ≥ (tr N)²/‖N‖_F². Row a contributes 1 + Σ_{b∈N(a)} ⟨v_a, v_b⟩². The v_b are orthogonal on non-free pairs, so they
form an orthonormal representation of the complement of F̄[N(a)], and Lovász's characterisation bounds the sum by θ_a. ∎

- **Polarity planes (B1, exact, every q).** In PG(2,q), y^⊥ ∩ x^⊥ is one point, so F[N(x)] is a matching and θ = 2. Every PSD
  decoder has rank ≥ V/3 = V/ω (dot and Hermitian forms, and the +1 charts).
- **Values** (`theta_local.py` recomputes three of them): P120 θ = 14, rank ≥ 8.000, factor ≤ 1.875 (tight); PG(2,4) dot θ = 2,
  rank ≥ 5.57 (global ϑ gives 6, so factor ≤ 0.889); PG(3,4) dot θ = 3, rank ≥ 16 = V/ω, factor 1.000. From the scouting run:
  PG(2,8), PG(2,16), PG(3,4) Hermitian and +1 charts all ≤ 1.000; PG(3,8) dot θ = 1 + 2√2, factor ≤ 1.207.
- **Erdős link.** A PSD factor is at most α(K, ω)/(ωK), where α(d, k) is the largest number of vectors in R^d with no k + 1
  pairwise non-orthogonal. A factor > 1 at ω ≤ 4 needs α(d, 4) > 4d. The search found α(4,4) ≥ 18, α(6,4) ≥ 27, α(7,4) ≥ 30 and
  α(8,4) ≥ 40 (possibly already known; compare Füredi–Stanley 1992).

## Q/F2 min-rank separation screen (expQA)

For field-line ports (θ = 1/n), the sandwich gives mr_Q/V ≤ γ_min ≤ mr_Q/(θV) = σ, and Theorem A gives rank ≥ θV in
characteristic 2. So the generic construction compresses iff σ < 1, and pays only if σ < 2/3. E8 P120 has σ = 0.60 at θ = 1/9.

| family (θ, h) | V | θV | rigorous mr ≥ (σ ≥) | exact decoder mr ≤ (σ ≤) | best numerics |
|---|---:|---:|---|---|---|
| PG(2,4) dot (⅓, 6), control | 16 | 5.33 | **7 (1.313)** C5-packing | **7**, so **mr = 7 exactly** | LM r = 7: 49/50; r = 6: 0/300 |
| PG(2,8) dot (⅓, 9) | 64 | 21.33 | 17 (0.797) C5 | 24 (1.125) | IRLS 25–28; seeded 24 |
| AG(2,8) +1 (⅓, 9) | 56 | 18.67 | 17 (0.911) C5 | 22 (1.179) | IRLS 24 |
| PG(3,4) dot (¼, 8) | 64 | 16 | 10 (0.625) α | 16 (1.000) | IRLS 16–20 |
| PG(3,4) herm (¼, 8) | 40 | 10 | 8 (0.800) C5 | 10 (1.000) | IRLS 10–13 |
| **PG(2,16) dot (⅓, 12)** | 256 | 85.33 | **52 (0.609) = α, Hoffman-tight** | 90 (1.055) | IRLS 98–103; seeded 91; LM r = 56 stalls at residual 90 = coverage law |
| PG(2,16) herm (⅓, 12) | 208 | 69.33 | 44 (0.635) | 72 (1.038) | seeded 72 |
| PG(3,8) dot (¼, 12) | 512 | 128 | 42 (0.328) | 141 (1.102) | seeded 141–142 |
| PG(2,32) dot (⅓, 15) | 1024 | 341.33 | 120 (0.352) | 354 (1.037) | IRLS 390 |
| E8 P120 (1/9, 9), calibration | 120 | 13.33 | 8 (0.600) | 8 (0.600) | rank 8 from random starts, 4/4 |

- **Upper bounds** are exact clique (frame) partitions N = Σ_C 1_C 1_Cᵀ; **lower bounds** hold over every field.
- **C5-packing bound.** mr ≥ Σ_i mr(G[U_i]) for vertex sets pairwise non-adjacent in the free graph, and **mr(C5) = 3 over every
  field** (rows 1 and 3 form a basis of a rank-2 N; the zeros N_41, N_53 force row 4 ∝ row 3 and row 5 ∝ row 1; then row 2 must
  vanish). This lifts PG(2,4) from α = 6 to 7, proving mr = 7, and the q = 8 families from 16 to 17.
- **α = 52 exactly for PG(2,16) dot.** The free graph is (255 vertices, degree 16, λ_min = −4) plus one isolated port, and
  Hoffman gives α ≤ 255·4/20 + 1 = 52, which CP-SAT and local search reach. So σ ≥ 0.609 and **no α-type argument can reach 2/3**
  there (Hoffman caps the α route at about 3/√q for PG(2,q)).
- **Negative results.** The Haemers (Hadamard) product bound is useless here (the least real M has rank 4 at q = 4 and 6 at q = 8,
  below α). **No perfect partitions:** PG(2,8) has none (CP-SAT optimal 24 > ⌈64/3⌉), and no perfect partition of PG(2,16) is
  invariant under the torus C_17.
- **Heuristic reason.** E8's separation is a Gram decoder whose required pairs form one residue class (even ⇒ 0). For F_q-lines
  with q ≥ 4 the required pairs carry q − 1 distinct residues, so at most one class can vanish: no Gram-type analogue exists.
- Only `verify.py` (the certificates) is replayed here; the IRLS, LM, CP-SAT and local searches are not.

## Search (h ≤ 17, ω ≤ 4 so that k ≥ 1/4)

| system | V | rank | factor V/(ωK) | label realisation |
|---|---|---|---|---|
| F4-18 (CP-SAT optimal in F4) | 18 | 4 | 9/8 | **none in PG(3,q) for any q**: 36/90 edges have a non-bipartite common free neighbourhood, but the polarity's is a matching (`obstruct.py`) |
| E6-27 (optimal in E6) | 27 | 6 | 9/8 | none in PG(3,4); none in PG(3,8)* |
| E8-40 (heuristic) | 40 | 8 | **5/4 (exact certificate)** | P120 restriction only (h = 9, u = 1, φ_eff = 0.56); none in PG(3,4); none in PG(3,8)* |
| BW16 (greedy) | 64 | 16 | 1.0 | — |
| Kneser coordinate labels KG(h, u) | — | — | ≤ 1 exactly (EKR) | native |

\* SAT homomorphism with port 0 fixed to each of 5 labels; complete only if Aut(T) is transitive. These SAT runs and the
heuristic subset searches are not replayed here.

## Cross-references

- `../../closures-2/field-line-large-ports/` (expNM): best exact decoder γ = 1 in all 8 polar families. Theorem A explains it
  in characteristic 2; Theorem B for PSD decoders.
- `../../closures-2/p120-rank8-rigidity/` (expQF): every rank-≤ 8 P120 decoder is Λ(G/2)Λ⁻¹. That is the E8 decoder that
  carries the 2-adic free cycles Theorem A2 requires.
- `../../closures-2/omega-h-floor-false/`: C ≥ ω·h is false. Theorem A is a rank bound uV/h, not a cost floor ω·h; it does not
  reintroduce the refuted screen.

## Next step and stop test (only if someone wants the real case closed)

1. Prove the real law at n ≥ 3 through A2: show that a polarity-graph decoder cannot carry a 2-adically negative free cycle
   without losing rank. **Stop:** a single LM/IRLS decoder below V/ω on any PG family refutes the real law.
2. The target needs α(d, 4) > 7d (about 1.4× better than E8-40) and an embedding into a polarity graph at h ≤ 17.
   **Stop:** K12, Λ16 or Leech subsets with ω ≤ 4 staying below 1.45·4·rank. Expected to stop.
3. **The only reopen route for large ports:** non-field labels with PR ∈ {1, …, u − 1}, whose free graph can be dense enough
   for an E8-type separation. Cheap first check: free-pair density ≥ 0.1; then expQA's E8-calibrated `irls.py` before anything else.
   The one rigorous gap left on field lines, σ ∈ [0.609, 0.667) for PG(2,16) dot, needs a non-α argument (e.g. a Schur complement
   on the Hoffman-tight coclique) and is not worth more than one session.

## Files and reproduction

Run each from this directory. Runtimes under a 4G cap, one worker.

| Script | Checks | Runtime | Output |
|---|---|---|---|
| `python3 -I -B trace_bound.py` | Theorem A's identity (N ⊗ J_u) ∘ G = blockdiag over F2 on PG(2,4) dot and P120, with random masked decoders; the 2-adic free triangle of G/2 | 0.4 s | `trace_bound.out`, `PASS_EXPQL_TRACE` |
| `uv run -q --no-project --python 3.12 --with numpy --with cvxpy python -I -B theta_local.py` | Theorem B values for P120, PG(2,4) dot, PG(3,4) dot (Clarabel SDP) | 2.1 s | `theta_local.out`, `PASS_EXPQL_THETA` |
| `uv run -q --no-project --python 3.12 --with numpy python -I -B certify.py` | exact certificate: F4-18, E6-27, E8-40 decoder ranks, ω = 4, factor 5/4; P120 labels reproduce E8-40's free graph, φ_eff = 0.5556 | 0.2 s | `certify.out`, `PASS_EXPQL_CERT` |
| `uv run -q --no-project --python 3.12 --with numpy --with networkx python -I -B obstruct.py` | F4-18 has no field-line realisation in PG(3,q) for any q | 0.6 s | `obstruct.out`, `PASS_EXPQL_OBSTRUCT` |
| `python3 -I -B price.py` | pricing through `kappa.py` (closures-2 copy, unchanged; it imports `../../framework-ceiling/`) | 0.1 s | `price.out`, `PASS_EXPQL_PRICE` |
| `qa/`: `uv run -q --no-project --python 3.12 --with numpy python -I -B verify.py` | expQA certificates: every partition decoder is an exact masked identity with the stated integer rank; every α / C5 bound re-checked; the σ-interval table | 0.4 s | `qa/verify.out`, `PASS_EXPQA_VERIFY` (also prints `ALL PASS`) |
| `lattices.py`, `fam.py`, `*.npy` | line systems; polar field-line families (closures-2 copy); the F4-18, E6-27 and E8-40 subsets | — | — |

Not included: the CP-SAT/heuristic subset searches (`maxsub.py`, `heur.py`), the LM fits (`lm.py`) and the SAT homomorphisms
(`hom.py`). They are search tools, and their results are quoted above.

Sources: expQL and expQA (with expNM and expQF as cross-references). `qa/` holds expQA's `verify.py`, `qa.py`, the certificate
JSON files, and copies of `fam.py` (expNM) and `minrank.py` (expNQ). Prepared with Claude Code (Claude Opus 5.5) assistance. Not
independently reviewed.

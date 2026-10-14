# P112 co-star nine-term masked decoders: none with 3–6 distinct normals

## Scope

**Exact non-existence result; no new κ, decoder or cost claim.** There is no nine-term masked
P112 co-star decoder whose terms use 3, 4, 5 or 6 distinct normals, over any field of characteristic
≠ 2, zeros and unequal left/right factors allowed. Normal sets with **k = 7, 8 or 9 distinct normals
are not covered** (k ≤ 2 cannot cover P112: no two co-stars do).

## Model and definitions

`engine/model.py` (independently rebuilt in `engine/verify_cert.py`, which pins PR #1320's line-order
SHA-256): the 120 E8 lines in #1320 order, Gram `G = v·w/4` (2 on the diagonal, 0/±1 off it); the frame
`F0 = {0,1,26,27,44,45,54,55}` is deleted, leaving **P112**. For a normal α (any of the 120 lines) the
**co-star** is `S_α = {r ∈ P112 : r ⊥ α or r = α}` (60 lines for live α, 56 for α ∈ F0).
A decoder is `M = Σ_{t=1..9} u_t v_tᵀ` on P112 × P112; term t carries a normal α(t) and
`supp u_t, supp v_t ⊆ S_α(t)`. The **capacity** `c_g` of a normal g is its number of terms
(`Σ c_g = 9`); group g contributes `P_g = D_g T_gᵀ` of rank ≤ `c_g`, with row/column supports `R_g, C_g`.
The **masked target** asks `M[r,s] = 2G(r,s)` on the 6,720 even ordered entries (`r = s`: 4;
`r ⊥ s`: 0); all other entries are free. Cases are (orbit of the normal set under the
172,032-element stabiliser H of F0, composition of 9 into k positive capacities); a set must be
**covering** (`∪ S_α = P112`, else some diagonal 4 is unreachable).

## Lemmas (all valid over any field, char ≠ 2; kills only set factor entries to 0)

*Ownership*: entry (r,s) is *alive* in g if `r ∈ R_g, s ∈ C_g`; it is *owned* by Γ if no group outside
Γ is alive there; it is *forced zero* if `r ⊥ s` or no group is alive there.

- **L0.** Every live r needs a group alive at (r,r), since `M[r,r] = 4 ≠ 0`.
- **L1 / L1Z.** For an orthogonal clique Q (L1Z: a *Z-clique*, all off-diagonal pairs forced zero)
  `M[Q,Q] = 4I`, rank |Q|; group g contributes rank ≤ `min(c_g, |Q∩R_g|, |Q∩C_g|)`, so this sum
  must be ≥ |Q|. "Initial L1" is this test with the untrimmed supports `R_g = C_g = S_g`, over all
  262,928 orthogonal cliques (the #1422 guard).
- **L2G / L2GZ (and L2, L2rank).** Let Γ be the groups alive on some row and some column of a
  (Z-)clique Q. If `Σ_Γ c_g = |Q|` then the stacked `D_Γ[Q], T_Γ[Q]` are square with product 4I,
  hence invertible. For r ∉ Q with `M[r,q]` forced zero and owned by Γ for all q ∈ Q,
  `D_Γ[r] T_Γ[Q]ᵀ = 0` gives `D_Γ[r] = 0`: delete r from every `R_g`, g ∈ Γ (columns symmetric).
  L2 is |Γ| = 1; L2rank: a clique of size `c_g + 1` owned by g alone needs rank `c_g + 1`.
- **L3 (generalises the #1320 crown lemma).** For `c_g = 2`, the *owned graph* has vertices x with
  (x,x) owned by g and edges x ⊥ y with (x,y),(y,x) owned by g. Then `D[x]·T[x] = 4`,
  `D[x]·T[y] = 0` in K². If y, z are both adjacent to x, `T[x]` annihilates `D[y]` and `D[z]`, so
  `D[y] ∥ D[z]`; adjacent vertices have independent D's. So the projective points of D two-colour each
  component: an odd closed walk is a contradiction (**L3odd**), and for x, y at odd distance `T[x], T[y]`
  are independent, so an owned zero row r with `M[r,x] = M[r,y] = 0` has `D[r] = 0` (kill).
  **L3GZ** (Γ of total capacity 2, forced-zero edges) is supported by the verifier but is used by no
  certificate; it was tried (engine5, not included) on the 23 survivors and closed none.

Every certificate is a list of kills plus a final L0/L1/L1Z/L2rank/L3odd contradiction, replayed with
plain Python sets by `verify_cert.py`, which imports neither `model.py` nor any engine.

## Closure table

| k | covering orbits | cases | initial L1 | propagation (rule at end) | Lemmas C+B′ |
|---|---:|---:|---:|---:|---:|
| 3 | 2 | 56 | 56 | 0 | 0 |
| 4 | 35 | 1,960 | 1,925 | 35 (L1 28, L0 7) | 0 |
| 5 | 567 | 39,690 | 39,140 | 550 (L1 345, L0 193, L3odd 11, L1Z 1) | 0 |
| 6 | 12,528 | 701,568 | 637,744 | 63,801 (L0 44,617, L1 19,029, L1Z 131, L3odd 24) | **23** |

k = 3, 4: `run_cert.py` + `engine2.py` (orbits from `cover.py`/`orbits.py`); k = 5, 6: orbits from
`group.py`/`canonfast.py`/`orbitsk2.py` (`orbk.pkl`), then `run_k.py`/`run_k6b.py` + `engine4.py`
(engine3 = L2G/L3, engine4 = + Z rules). Propagation kills used: k=5 L2G 101,808, L3 786.

## The last 23 cases (k = 6, all capacities (2,2,2,1,1,1)) — `lemmas23/`

Each has a pairwise non-orthogonal A2 triple {a,b,c} of capacity 2 and three mutually orthogonal
capacity-1 normals {d,e,f} orthogonal to the triple. `class(a)` = live lines ⊥ a but not ⊥ b, c.
No `S_b, S_c` meets class(a), so on class(a) only `P_a` (rank ≤ 2) and the rank-1 terms x ∈ {d,e,f}
(supported on `W_x = {r : r ⊥ x}`) act; on a cross-class orthogonal pair only rank-1 terms act.

- **Lemma C (SAT, `check_C.py`).** Booleans for the supports of `u_x, v_x`. Clauses: on each
  cross-class orthogonal pair the number of nonzero products `u_x[r] v_x[s]` is never exactly one;
  for each orthogonal triangle T inside a class, `det(4I − Q_T) = 0` (P_a has rank ≤ 2) while
  `det = 64 ≠ 0` unless some loop/2-/3-cycle of `Q = Σ u_x v_xᵀ` has all entries nonzero. Claim:
  not every class is served (u_x, v_x both nonzero on it) by ≥ 2 terms. UNSAT in all 23 cases with
  CaDiCaL 1.5.3 and Glucose 4.
- **Lemma B′ (Nullstellensatz, `check_B.py`).** In class(a), `P_a` plus at most one rank-1 term
  `u vᵀ` on `W_x` cannot give 4 on the diagonal and 0 on orthogonal pairs. Gauge: two orthogonal lines
  r0, r1 ∉ W_x see only P_a, so `U_r0, U_r1` are independent; GL₂ fixes them to e1, e2; rank-1 scaling
  splits `u[w0] ∈ {0, 1}`. 23 × 3 classes × 3 terms × 2 = **414** certificates `Σ g_i f_i = 1`, cofactor
  denominators all powers of 2, so infeasible over any field of characteristic ≠ 2.

Lemma C leaves a class served by ≤ 1 term; Lemma B′ excludes that block. All 23 cases are infeasible.

## Addendum: every normal set containing F0 is closed (`f0/`)

These are the k = 8 set F0 with caps (2,1⁷) for any choice of cap-2 normal, and the k = 9 sets F0 ∪ {x} with caps 1⁹ for all 112 live x. Support-level rules cannot close them; the following closes them, valid over any field of characteristic ≠ 2.

- **Structure.** The 112 live lines split into 14 blocks of 8. Each block is two orthogonal 4-cliques whose lines share the same 4 F0 seers, giving 28 cliques in all.
- **Triangle lemma.** Take three cliques that pairwise share exactly one F0 pair P, together with matched 2-subsets forming an orthogonal 6-set. Assume two of the cliques are "pure-like" and the third has at most 3 other terms. Then the system is infeasible.
  - The idempotent E = ¼V[A]ᵀU[A] and its complement force the third clique's projection Π to be 0.
  - Then 4I − (nilpotent) has rank 4 > 3, a contradiction.
  - Certificate: a Nullstellensatz identity Σgᵢfᵢ = 1 (`triangle_cert.json`; 48 variables, 57 equations, 849 cofactor terms, degree ≤ 6, all denominators powers of 2).
- **Case A (F0, cap 2 on g0).** For every g0, some triangle has two pure cliques and one g0-clique.
- **Case B (F0 ∪ {x}).** A SAT encoding on the support of the x-term is UNSAT for all 112 x under both CaDiCaL and Glucose. It uses two clause families:
  - FR: completely orthogonal clique pairs seen only by x;
  - TRI: one clause per triangle and choice of third clique.
  - Controls: dropping either clause family makes the encoding SAT (`negctl.py`).

Run `python3 f0/check_F0.py`. It needs python-sat, rebuilds the E8 lines, imports nothing else from the package, and takes about 3 s. `f0/triangle_cert.py` regenerates the certificate with passagemath-singular.

Other k = 7, 8, 9 normal sets remain open (see Discussion #1267).

## Not covered

- k = 7, 8, 9 distinct normals (not attempted here).
- Decoders whose factors are not confined to the co-star of their normal (non-co-star supports).
- Masked decoders of rank ≤ 8 (fewer terms): handled by #1212 (rigidity), not here.

## Reproduce

Python ≥ 3.12 + numpy; cap memory (e.g. `systemd-run --user --scope -p MemoryMax=6G`).
Times: 8-core Ryzen AI 7 PRO 350, Python 3.14 + numpy 2.5; `lemmas23` under `uv` Python 3.12.

```sh
cd engine
python3 verify_cert.py cert3.pkl cert4.pkl        # stdlib only, 94 s
#   cert3.pkl : k = 3 reps 2 certificates verified 56 / 56
#      covering 3-sets: 1120 ; orbits 2 ; reps given 2 ; sets not in a rep orbit: 0
#   cert4.pkl : k = 4 reps 35 certificates verified 1960 / 1960
#      covering 4-sets: 245224 ; orbits 35 ; reps given 35 ; sets not in a rep orbit: 0
python3 verify_k.py orbk.pkl certk5.pkl            # 4 workers, 128 s
#   k=2: 8 orbit reps; sizes sum to C(120,2); all extensions covered (6s)   [k=3: 40, k=4: 296 likewise]
#   k=5: 2701 orbit reps; sizes sum to C(120,5); all extensions covered (128s)
#   certk5.pkl: k=5 covering reps 567, (rep,composition) pairs 39690, certificates verified 39690, open 0
cd ../lemmas23
uv run -q --no-project --python 3.12 --with python-sat --with sympy python check_C.py    # 1 s
#   OPEN list matches ../engine/k6_open.pkl (23 cases)
#   23 lines '<normals> <caps> class lines 75|78 cross pairs (…) cadical UNSAT glucose UNSAT'
#   LEMMA C (each class served by >=2 terms is impossible): PROVED for all 23
uv run -q --no-project --python 3.12 --with python-sat --with sympy python check_B.py    # 16 s
#   LEMMA B' certificates verified: 414 (all cofactor denominators are powers of 2) (23 cases x 3 classes x 3 terms x 2 scaling cases)
```

**k = 6**: `engine/regen_k6.sh` (`k6chunks/` absent) rebuilds `cliques.pkl` (byte-identical), runs
filter + engine4, writes 40 chunks of 701,545 replay jobs, checks `k6_open.pkl` unchanged, replays all.
Original run, 15 workers: 165 s filter + 1,610 s engine, **~30–45 min** end to end (`JOBS=4`: ~1–2 h).
Expected: `63824/63824 … {'L0': 44617, 'L1': 19029, 'L1Z': 131, 'L3odd': 24, 'SURVIVE': 23}`, then
`k=6 certificates replayed: 701545`; original chunk listing `sha256sum c??.pkl | sha256sum` =
`a597fab2c9b1f8a9ae1e3fe914a9b3d8c3b1ab41ea1aa3d9f0d55cf29e26f22a` (the engine's 20 s Z-scan limit
means a much slower machine could change bytes; the counts are the invariant). Here: the original
chunks replayed with the packaged `verify_cert.py` (701,545 pass, 16 s; step-stripped certificates
are rejected); engine4 re-closed 60 sampled cases identically and left the 23 open. A full `JOBS=4`
rerun was stopped at 32,000/63,824 jobs (3,338 s, heavy load from other jobs) to respect a 1-hour
budget; its 16 checkpoints matched the original log exactly. k = 6 orbit completeness
(`ORBIT_ONLY=1 python3 verify_k.py orbk.pkl certk6.pkl`, 2,606 s originally) was not repeated.

SHA-256: `cert3.pkl` ed8b38bb…d706, `cert4.pkl` be38d5f3…7e197e, `certk5.pkl` 1cb80c84…2c4d,
`orbk.pkl` 774b1d36…0e49, `k6_open.pkl` 14f2d1ec…f5f6.

Not committed (size): `certk6.pkl` (~490 MB), `certk6_filter.pkl`, `k6chunks/` (40 × ~12 MB),
`cliques.pkl`, `group_elems.npy` — all regenerated by the scripts above. Producers needing Sage
(`gbB.py`, `gbB2.py`, `liftB.py`) or pycryptosat (`sat1.py`, `sat4.py`) are included for provenance.
Worker counts are `JOBS` (default 4; the original runs used 15). `verify_k.py` forces the `fork` start method.

## Credits

Builds on Hanhong Zhao / @Hope7Happiness (#1267; #1320 crown lemma, generalised here as L3),
@jeremydiviney (#1422 clique demands, the L1 guard), @burntout76 (#1212 rigidity),
@DreamingOfClouds (#717) and @themilkman216 (#968).

Prepared by Philippe Van Velzen with substantial Claude (Anthropic) assistance; parallel sub-agents did the computations, the main session replayed the k≤4 certificates and both 23-case checkers independently. Apache-2.0.

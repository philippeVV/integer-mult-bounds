# Conjecture Θ, partial proof (expXR): second order on the whole flat manifold, an exact global class, and a refuted fallback

Status: theorems plus numerical evidence, author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`.
Multipliers are over κ = 0.0017. Builds on `../../limits-1/complex-dual/` (expXK Lemma 1, Thm 2, Prop 4/5; expXG).

## Verdict

| Item | Statement | Status |
|---|---|---|
| Question | Conjecture Θ (girth 8, `../../limits-1/complex-dual/`): every MDS half-rate complex Tanner code at girth ≥ 8 has a gauge with Θ̂₃ ≤ 1/8 | **open**; this lane proves it at second order everywhere on the flat manifold and exactly on a global class |
| **Theorem L** | at every two-sided-flat code C₀ (any even degrees, girth ≥ 8), the second-order code/gauge max–min model is ≤ 0. In particular Γ ≻ 0 ⇒ the normal Schur complement S ⪯ 0. Mechanism: the inertia count **n₊(K) ≤ E − 1** | **proved** (Lemmas A, B); checked L1–L4 |
| **Lemma A** | for symmetric contractions X, Y the 2E×2E form H₀ = [[4Y, 3(I+YX)], [3(I+XY), 4X]] has n₊(H₀) ≤ E; it is nonsingular with inertia (E, E) when X ≻ −I or Y ≻ −I, and singular iff XY has eigenvalue e^{±2iθ}, cos 2θ = −1/9 (sharp) | **proved**; L2a–c |
| **Lemma B** | H₀(1,1) = 10(1,1) and (1,1) ⟂ every realisable leverage deviation, so n₊(K) ≤ E − 1 = the number of non-trivial gauge directions | **proved**; L3 (W3, W4M, W5; random, real, Fourier codes) |
| Boundary Γ ⪰ 0 | the refutation route "flat code with Γ ⪰ 0 and a positive normal Schur direction" is impossible for Γ ≻ 0, and for singular Γ ⪰ 0 whenever H₀ is nonsingular | **proved** (§2.4 of the source) + adversary: best re-evaluated values 0.12498403, 0.12498779 < 1/8 |
| **Global reduction** | every MDS code is (C₀ flat, δ) by per-vertex Barthe flattening, so Θ ⇔ Φ_{C₀,δ}(1) ≤ 1/8 for all flat C₀ and δ, where Φ(t) = inf_x Θ̂₃(C₀; points e^x, lines e^{x+tδ}) | **proved** |
| **Conjecture M** | Φ is non-increasing on [0, ∞). M ⇒ Θ; Theorem L is exactly M at t → 0⁺ | **open**; stop test negative (evidence below) |
| **Theorem P** | at the point-flat gauge, with ε = line leverage − ½: Θ̂₃ ≤ 1/8 ⇔ Σ_p [¼Σ_{i<j} X^p_{ij}ε_iε_j + 6Σ_{i<j<k} Re(Π^p_{ij}Π^p_{jk}Π^p_{ki})ε_iε_jε_k] ≤ (εᵀXε)²/(2E). Exact, any degrees | **proved**; L8 (moment error ≤ 6e-17) |
| **P1** | if after point flattening no point lies on two non-flat lines, Θ̂₃ = 1/8 exactly, so Θ̂* ≤ 1/8. Line codes and point codes otherwise arbitrary. Strictly contains the two-sided-flat class (W3: a spread of 10 arbitrary lines; W5: ≥ 19) | **proved** |
| **P2** | at most two non-flat lines per point with ε_iε_j ≤ 0 at each such point also gives ≤ 1/8. Mirror statements (lines flat) hold too | **proved** |
| **Linear-gauge fallback** | "x = −Γ⁻¹Cδ keeps Θ̂₃ ≤ 1/8 on Γ ≻ 0 bases" | **false**: W3, Dirichlet(4) codes, λ_min Γ = 5.8e-3, ‖δ‖_rms = 0.6 gives **0.1253568** (dense, independent evaluation) while the true inf over gauges is 0.1227497 (L5). On SIC it fails from ‖δ‖_rms ≈ 3 (0.1295) |
| One-sided gauges | point-flat and line-flat gauges each exceed 1/8 on the W3 Γ-boundary adversary codes and on SIC+noise (0.12505–0.12513); a proof of Θ needs interior gauges | numerical (`pl_test.py`) |
| Multiplier | **×0** (ceiling work; no κ is produced) | |

## Setting

Girth ≥ 8, even degrees, MDS half-rate complex local codes. By Barthe/Forster, A_v = e^{f_v}A_v⁰ with A_v⁰ flat (leverages ½),
unique up to phases. With u = f_P, v = f_L, the code is gauge-equivalent to C₀ = (A_v⁰) with points scaled by e^x and lines by
e^{x+δ}, δ = v − u modulo vertex constants (E − V + 1 directions, the 81 "normal" directions expXK found at SIC on W3). expXK
Thm 2 gives Θ̂₃(C₀, x = 0) = 1/8 for every flat C₀. Degree 4: flat [4,2] codes ↔ w in the simplex (Bloch vectors), so degree-4
codes ↔ (w ∈ Δ^V, δ ∈ ℝ^{E−V+1}) (`model.py: flat_code`).

## Theorem L (proof outline)

- Second-order model at a flat code: α = ½(I−X)u, β = ½(I−Y)v with X = 4|Π_P|^{∘2} − I, Y = 4|Π_L|^{∘2} − I (block diagonal,
  symmetric, zero diagonal, row sums 1, ‖X‖, ‖Y‖ ≤ 1 by Schur's product bound). Then f(u, v) − 1/8 = (4/E)·K(u, v) + O(3) with
  K = (1/32)(α,β)ᵀH₀(α,β). First-order terms vanish identically, so moving along the flat manifold is free at second order and the
  max–min is sup_δ inf_x K(x, x+δ). The diagonal block K(x, x) is expXK's Γ (L1, L1b; model/exact ratio 0.9998–1.0017).
- **Lemma A.** ker H₀ is invariant under T(α,β) = (Xβ, Yα); a kernel eigenvector forces 3λ² + 4λ + 3 = 0, λ unimodular and
  non-real with Re λ² = −1/9; contractivity rules this out when X ≻ −I (or Y ≻ −I). H₀ is nonsingular on that convex set ∋ (0,0),
  where its inertia is (E, E); limits only drop n₊, so n₊ ≤ E for all contractions. Sharpness: reflections at the critical angle.
- **Lemma B.** K = (1/128)ZᵀH₀Z with Z = diag(I−X, I−Y); c = (1,1) has H₀c = 10c and c ⟂ range Z, so any H₀-positive subspace
  of range Z extends by c: dim P + 1 ≤ E. (Raising every leverage is positive but not realisable: Σ_{e∈v} a_e = deg(v)/2.)
- **Theorem L.** If Γ ≻ 0 on 1^⊥ then n₊(Γ) = E − 1 and Haynsworth gives n₊(S) = n₊(K) − n₊(Γ) ≤ 0, so
  Θ̂*(C₀ ⊕ tδ) ≤ 1/8 + (4/E)t²δᵀSδ + O(t³). If Γ has a negative direction, Θ̂* < 1/8 nearby (expXK Prop 5). L4 finds exactly
  E − V + 1 negative Schur eigenvalues at SIC and generic flat codes on W3 and W5. ∎

## Theorem P (proof outline)

All point codes flat (a ≡ ½, |Π_p|^{∘2} = ¼(I + X_p)); line codes arbitrary with b = ½ + ε, Σ_{e∈l}ε_e = 0. Evaluating each
term of expXK Lemma 1: line traces take their flat values whatever the line code; point traces pick up εᵀ|Π|^{∘2}ε and the
cubic T_p(ε) = Σ Π_ijΠ_jkΠ_ki ε_iε_jε_k; all Σε terms vanish. Result:
E m₂ = 3E/16 + ¼q, E m₃ = 5E/32 + (3/8)q + C, with q = εᵀXε and C = Σ_p 6Σ_{i<j<k} Re(Π_ijΠ_jkΠ_ki)ε_iε_jε_k (repeated-index
terms cancel). The Markov bound at s = ½ + 2q/E (< 1 since |q| < E/4) gives Θ̂₃ − 1/8 ≤ [q/8 + C − q²/(2E)]/E·(1−s)⁻². ∎
P1: ≤ 1 non-zero ε per point ⇒ q = C = 0. P2: ≤ 2 per point ⇒ C = 0 and q ≤ 0 when ε_iε_j ≤ 0, since X_ij = 4|Π_ij|² ≥ 0. ∎

## Conjecture M: evidence and stop test

- **Differential form.** Danskin plus gauge stationarity give tΦ′(t) = (d/ds) f(su, sv)|_{s=1} at the optimal gauge: M says that
  pulling every vertex code radially towards its own flat representative never decreases Θ̂₃ at a gauge optimum. Confirmed
  to 1e-6 (L6: FD −1.91456402e-02 vs ⟨g_v, δ⟩ = −1.91456403e-02).
- **Rays.** 39 rays (W3 SIC; Γ ≻ 0 and Γ-indefinite bases at Dirichlet 0.1/0.3/2/4; W5 generic flat; the W3 Γ-boundary adversary
  code), t ∈ [0, 6]: all monotone, none above 1/8. The ray through expXK's W5 adversary best code (Θ̂* = 0.12499917) is
  monotone too (0.12499995 at t = 0.25, 0.12499917 at t = 1). L7: a precisely re-solved ray through a base with
  λ_min(Γ) = 5.2e-4 is monotone (0.1248402 → 0.1055097).
- **Stop test (adversary on J(δ) = Φ(δ) − Φ(0.8δ), `m_adv.py`).** 9 runs × 30 ascent steps on W3 (SIC, Dirichlet(2), 97% of the
  way to ∂{Γ ⪰ 0}; ‖δ‖_rms = 0.05 … 2): max J ∈ [−6.8e-3, −9.5e-5], every re-solved ray monotone. **Negative.** 11 queued runs
  were stopped for budget.
- **Caveat.** The inner problem has several local minima for ‖δ‖ ≳ 0.3, so every computed Φ is an upper bound on the true Φ.
  A monotone sub-1/8 upper-bound curve is evidence, not proof.
- **No per-vertex identity.** The vertex terms of the second-order radial derivative are all ≤ 0 at SIC (0/1600 positive) but
  16/1600 are positive on Dirichlet(4) bases (`pervertex_q.py`), so any proof of M has to be non-local (through the gauge
  multiplier).

## Next step and stop test

1. **Interior gauges.** Theorem P and its mirror cover the two endpoint gauges; the endpoints alone fail. Target the family
   x_s = (1−s)·(point-flat) + s·(line-flat): first the expansion of Θ̂₃(x_s) in s at s = 0 (needs the line-side Lemma 1 terms at
   non-flat lines). **Stop:** a code with min_s Θ̂₃(x_s) > 1/8 while Θ̂* ≤ 1/8 (linear interpolation not enough either).
2. **Conjecture M.** **Stop:** a precisely re-solved ray with Φ′ > 0 (that would also be the first candidate counterexample to Θ).
3. The linear-rule fallback is closed: it certifies only on a δ-ball whose radius shrinks with λ_min(Γ) (r ≈ 0.6 at
   λ_min = 5.8e-3; r ≈ 2–3 at SIC).

## Files and reproduction

All scripts are in `xr/` and import expXK's `xk_sparse.py`, `gauge_hess.py` and expXG's `xg_core.py` from
`../../limits-1/complex-dual/{xk,xg}/` (graphs W3, W4M, W5 from the `graphs_*.json` there); `model.py` sets that path.

| Script | Checks | Runtime (one worker) | Output |
|---|---|---|---|
| `uv run -q --no-project --python 3.12 --with numpy --with scipy python -I -B check_xr.py` (slow) | L1–L4, 19 checks: second-order model and Γ; Lemma A (roots, sharpness, inertia); Lemma B on 9 codes; Theorem L (Schur ⪯ 0, E − V + 1 negatives) | 97 s | `check_xr.log`, `PASS_EXPXR` |
| `uv run … --with numpy --with scipy python -I -B check_xr2.py` (slow) | L5 linear-rule counterexample (dense, independent), L6 Danskin, L7 precise ray, L8 Theorem P on W3/W5 incl. the P1 class | 51 s (about 5 min under load in scouting) | `check_xr2.log`, `PASS_EXPXR2` |

Both replays reproduce the scouting logs line for line (apart from the timing). Included but not replayed (scouting sources
with their logs): `pointflat.py`, `pl_test.py` (one-sided gauges); `lin_adv.py`, `lin_sphere.py`, `lin_*.log`,
`linsph_W3_4_5_r0.6.npy` (input of L5); `m_adv.py`, `batch_m.sh`, `batch_m[AB].log`, `madv_*.log` (Conjecture M adversary);
`pervertex_q.py`; `inertia*.py`, `h0test*.py` (inertia scans, sharpness search); `mk_bdry.py`, `adv_bdry_s*.log`,
`eval_s1*.log` (Γ-boundary adversary, run through `../../limits-1/complex-dual/xk/outer_xk.py` and `evalbest.py`); `ray.py`,
`ray2.py`, `ray_code.py`, `ray2_*.log`, `raycode_W5.log` (ray scans). Paths inside these sources still point to the scouting
layout (`../expO/run`, `../expXK/`).

Source: expXR. Prepared with Claude Code (Claude Opus 5.5) assistance. Theorems L, P (with P1, P2), Lemmas A, B and the global
reduction are proved here; Conjecture M and the adversary results are numerical. Not independently reviewed.

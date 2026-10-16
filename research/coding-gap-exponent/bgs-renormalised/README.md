# BH⊗BGS with a renormalised starting dual

**BH⊗BGS: the guaranteed exponent is ≈ 0.0056 (all N), not 0.0020. This is a renormalisation of BGS's own starting dual.**
- BGS App. C/F start the amplification from the α-extension A* of the seed. For the Braverman–He seed (735 edges, 157 sessions, all at distance 5), A* is the 735 edges plus 314 unit-capacity pendant leaves.
- The dual D_α they use keeps the core weights 1/785 and gives the leaves weight 0.077× an edge weight. It has slack Σ dist = 1.031 and is never rescaled. The result is γ₀ = 1+ε/2 = 152/147 and c = l/w = 66.94, so all-N η = ln γ₀/(4 ln(c+1)) = 0.00198.
- On the same graph, the uniform dual w ≡ 1/1099 is feasible with Σ dist = 1 exactly. It gives z = 1049/1099, γ₀ = 1099/1049 and c = 7, so **η = ln(1099/1049)/(4 ln 8) = 0.00560** for all N, and 0.0112 along BGS's sequence. With the printed Thm 1.1 constant ln γ₀/ln(16c₁²), the sequence value is 0.0067.
- **Optimality on this graph.** Any positive dual with l ≤ c·w has z ≥ a + 2(1−a)/c with a = 147/157, by weak duality against the core flow, and c ≥ 7. So uniform is the best dual for every rounding level. LP: HiGHS cutting planes, which agree to 10 digits.
- **Why the honest optimal dual does not help.** The core MCF is exactly 147/157 (rational flow certificate), and every optimal dual of the leaf graph puts weight 0 on all leaves. That is BGS's own "weight 0" obstruction, so a non-optimal dual is unavoidable.
- **Hypotheses.** Equal demands, distinct terminals, positive weights, feasibility and c_m ≥ c_k, c_n all hold. The coding rate stays ≥ 1 because each leaf carries one F₉ symbol per slot. App. F's proof uses D* only through standard form, NC/z(D*), l and w (Thm 3.6, Lemma 3.9, w_i = w_{i−1}²). The printed Thm 3.10/1.1 cannot be applied directly, since Lemma 2.2 would contract the leaves. **Cite as "App. F's proof with D* = uniform".**
- **Further gain under the same standard form.** Keep one terminal at each of the 133 terminal vertices and leaf only the other 181 terminals. This gives γ₀ = 483/458 at c = 7, so η = **0.00639** (all N) and 0.0128 (sequence).
- **Consequence for the repo.** Upstream's κ_sel = 0.00167 is about 0.26–0.30× of the guaranteed undirected coding-gap exponent. It is not comparable to it.
- Scripts: exact Fraction checker plus HiGHS LPs (expXI). Not independently reviewed.

Prepared with Claude Code (Claude Opus 5.5) assistance; not independently reviewed.

Run: `python3 cert.py` (exact Fractions; prints the three dual variants). LP scripts (`task1.py`, `tradeoff.py`, `cp.py`) need `highspy`.

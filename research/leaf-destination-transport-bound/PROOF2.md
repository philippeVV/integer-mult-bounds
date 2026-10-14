# Exact all-size value: E(n) = max(1, floor(log2 n))

The domain, notation, Facts F1–F2 and Lemmas S, L2, L1 are as in `PROOF.md`. A call of size m
has B = its first b = floor(m/2) factors and A = its last a = ceil(m/2) factors. Its scratch
levels are numbered 1 (z_K), 2 (the children's z), and so on. fl(x) = floor(log2 x).

For a call K, e(K) is the maximum mass of any register at any cut b >= a inside K, where cut a
follows an injection at a leaf of K. Put e(m) = max over calls of size m, both signs. Then
E(n) = e(n). We use fl(floor(k/2)) = fl(k) - 1 and fl(ceil(k/2)) = fl(k+1) - 1 for k >= 2.

**Theorem 1.** e(m) <= max(1, fl(m)) for all m.
**Theorem 2.** For every n >= 2 there is an explicit interval of mass fl(n).
So E(n) = max(1, fl(n)) for all n. In particular E(2^h) = h.

## 1. Two new lemmas

**Lemma I (inverse).** word(k, -s) is word(k, s) reversed, with every leaf sign negated. Since
a leaf has x != y, the leaf (x, y, f, -s) is the inverse of (x, y, f, s). So the sign -s call is
the operator inverse of the sign s call, leaf by leaf.
*Proof.* Induction. The sign - call reverses the phase list and negates each phase's sign. ∎

**Corollary I.** Let C and C' = C^{-1} be calls on the same registers. Let column e come from an
injection at a leaf of C, followed by C's suffix. Then any prefix of C' maps e to a forward
transport inside C or inside C'. That transport starts at an injected leaf destination. So its
mass is at most e(|C|).
*Proof.* Let C have leaves L_1..L_N. A prefix of C' of length j is (L_{t+1} ... L_N)^{-1} with
t = N - j. If t >= i, the product is the segment L_{i+1}..L_t of C. If t < i, the inverse
L_i^{-1} does nothing to e_{y_i}, because it reads x_i. The rest, L_{i-1}^{-1} .. L_{t+1}^{-1},
is a segment of C' that starts right after its leaf L_i^{-1}, whose destination is y_i. ∎

**Lemma X (exit support).** Inject at a leaf of a call K of size k. At K's end, let S be the
support of the column and S_z its scratch part. Then |S| <= fl(k) + 1 and |S_z| <= fl(k).
*Proof.* Induction; k = 1 gives S = {y}. Let φ be the phase of K that contains the leaf. The
L1 case analysis gives the following. Phases after φ never change the entries on scratch levels
>= 2, because they are complete calls and F2 applies. They also leave z_K unchanged. They can
only set y_K, and only when φ = U' (sign +) or φ = U (sign -), where y_K = -A c_z.
- φ of V-type (size a): S = S(φ) and S_z = S_z(φ). So |S| <= fl(a) + 1 <= fl(k) + 1 and
  |S_z| <= fl(a) <= fl(k).
- φ of U-type (size b): y_φ = z_K, so S ⊆ S(φ) ∪ {y_K} and S_z = S(φ). So
  |S| <= fl(b) + 2 = fl(k) + 1 and |S_z| <= fl(b) + 1 = fl(k). ∎

Both bounds are attained. Call the U-phase that keeps y "y-keeping": U' for sign +, U for
sign -. Follow the y-keeping U-phase at every level (these alternate in sign). The resulting
leaf gives a support of all fl(k) scratch levels plus y. Each step prepends -A_i, so the entries
are ± monomials of strictly increasing length, from the unit at the deepest level up to y.

## 2. Proof of Theorem 1

Use strong induction on m. If both cuts lie in one child, use e(child). Otherwise K (size m >= 2)
is their lowest common call. Cut a lies in phase φ_α and cut b in a later phase φ_β, and
M = Pre_β(b) ∘ T_mid ∘ Suf_α(a). Let q be the column at φ_β's start. All its entries have mass at
most 1 (L1, F2). Registers outside φ_β's footprint keep that bound. By L2, a footprint row has
mass |Σ_v Pre[w,v] q_v| <= #{v in footprint : q_v != 0} (*).

| K | (α, β) | argument | bound |
|---|---|---|---|
| + | (U,V), (U',V') | T_mid = I. (*) with q = exit(φ_α), Lemma X | fl(b)+1 = fl(m) |
| − | (U,V') | same | fl(m) |
| + | (U,V') | Row y gets (A + Pre[y,z]) c_z, which lies in {A c_z, 0} by S. Other rows have Pre[w,y] = 0, so the count is at most \|exit(U)\| | fl(m) |
| ± | (U,U') | T_mid touches only y_K, which U' does not read. Corollary I with C = U | max(1, fl(b)) |
| ± | (V,V') | T_mid is U' (+) or U (−) and reads x_K = 0, so it acts trivially. Corollary I with C = V | max(1, fl(a)) |
| + / − | (V,U') / (V,U), (V,U'), (V',U') | x_K = z_K = 0 (T_mid reads only zero entries). (*) counts S_z(exit V-type) ∩ levels(U) | fl(a) |

Here fl(b) + 1 = fl(m), fl(b) = fl(m) - 1 and fl(a) = fl(m+1) - 1 <= fl(m). Every case is at most
max(1, fl(m)), and fl(m) <= fl(n). ∎

**Exact case values.** For all m <= 80 and both signs, the census matches the upper bounds in the
last column exactly:
- fl(m) for (U,V), (U,V'), (U',V') (sign +) and (U,V') (sign −);
- max(1, fl(m) - 1) for (U,U');
- max(1, fl(m+1) - 1) for every case that starts in a V-phase.

So the maximum is attained in the four "source-dirty" cases for every m. It is also attained in
the V-start cases exactly when m + 1 is a power of two. The UU' case never attains it once m >= 4.
The lower bound is proved below for (U,V) with sign +. The other three source-dirty cases are
census-only. Lemma P below fails for sign − (see §3).

The old loss of 1 came from counting ceil(log2 |φ_β|) scratch levels. Lemma X instead counts the
dirty registers that φ_α can export: at most fl(b) + 1, with b the floor half. Only B-phases
export to a footprint that holds the extra "+1" source register.

## 3. Proof of Theorem 2

**Lemma P.** Take a call of sign + and size k. Let t be the cut after its V-chain leaf, where we
follow V at every level; this is the first leaf of factor f_{k-1}. Then Pre(t)[y, x] = P_K and
Pre(t)[y, z_ℓ] = f_{k-1} ... f_{k - ceil(k/2^ℓ)} for every scratch level ℓ = 1..ceil(log2 k).
*Proof.* Induction. Write K = U V U' V' and note t lies in V. Then T = U is z += B x, and y is
untouched. So Pre_K[y, x] = Pre_V[y, z] B and Pre_K[y, z] = Pre_V[y, x_V], since z = x_V. The
deeper levels are V's scratch levels, which the induction covers. ∎
(This fails for sign −. There T contains V, which adds A z to y, and V' cancels it in row y.
The check shows the failure.)

**Family.** Take n >= 2, the top sign +, and H = fl(n). Inject at the leaf of factor f_0 reached
inside U along the y-keeping chain U, U', U, U', .... Read row y at the cut of Lemma P inside V.
By Lemma X's tight case, q = exit(U) is supported on z_K = level 1 and on levels 2..H. Its
entry at level ℓ is (-1)^{H-ℓ} f_{floor(n/2^ℓ)-1} ... f_1. Lemma P applies to V, with the levels
shifted by 1, and H - 1 = fl(b) <= ceil(log2 a). It gives

  y = Σ_{ℓ=1}^{H} (-1)^{H-ℓ} · f_{n-1} ... f_{n-ceil(n/2^ℓ)} · f_{floor(n/2^ℓ)-1} ... f_1.

The original y entry is 0. The H words have pairwise distinct right factors over the B-letters,
because their lengths differ. So nothing cancels, and the mass is exactly H = fl(n). ∎

For n = 2^h this gives Σ_{j=0}^{h-1} (-1)^j (f_{n-1} ... f_{n-2^j})(f_{2^j-1} ... f_1). For
example, n = 16 gives f15 - f15f14f1 + f15..f12f3f2f1 - f15..f8f7..f1.

## 4. Checks (exact free algebra)

- `check2.py ok LO HI`, run as 1–56, 57–64, 65–74 and 75–80, checks all of the following for each n, with both orientations and every
  injection-cut pair:
  - INV (Lemma I);
  - L1 (monomial exit, zero source);
  - Lemma X (≤ and tight);
  - the case table equal to the closed forms;
  - E(n) = max(1, fl(n));
  - Lemma P (sign +);
  - the family mass fl(n) with distinct exit lengths.
- The family check also runs alone for n <= 512.
- `fam_closed.py 2 512` checks the family's row y against the displayed formula exactly.

Mutation controls (`controls2.log`):
- unreversed inverse: INV, L1, CASE and E fail; E(7) = 54.
- swapped split (B = ceil half): X, CASE, P, FAM and E fail; E(6) = 3 > fl(6). This shows that
  the floor-half export count matters.
- dropped top restoration: INV fails.
- The swapped-split family misses the closed form.

Logs: `check_ok_*.log`, `fam_closed.log`, `controls2.log`.

Not covered: anything outside PROOF.md's scope, namely internal leaf intermediates, multiple or
accumulated errors, rounding, aliasing and kappa.

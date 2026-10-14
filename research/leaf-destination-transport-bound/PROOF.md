# All-size bound: destination-supported transport is at most h = ceil(log2 n)

Domain: the one in `research/leaf-destination-transport`. PR628 `compose` (mid = n//2, B = first
half, A = second half), with sign + phases U+,V+,U-,V- and sign - phases V+,U+,V-,U-. Here
U = (x -> z, B), V = (z -> y, A) and z is the depth scratch. Factors are free noncommuting
contractions. The mass |p| is the sum of absolute coefficients, and it bounds the norm of that
entry. A unit error goes into the destination of a completed leaf at cut a. E(n) is the maximum
mass of any register at any cut b >= a, over both orientations.

**Theorem.** E(1) = 1, and E(n) <= ceil(log2 n) for n >= 2.

This proves max(1, floor(log2 n)) for powers of two. Combined with PR1416's lower family
(h - 1), it leaves the maximum in {h - 1, h}.

**Facts.** (F1) Inside a call K, the source x_K is never written and the destination y_K is never
read. (F2) A complete call K of sign s acts on any column exactly as y_K += s P_K x_K, where
P_K is the ordered product of its factors. Its scratch is restored exactly.

**Lemma S and Lemma L2, proved together by induction on k.** Let Pre_K(t) be the map from K's
start to any cut t inside K.
- (S) Pre_K(t)[y, x] is 0 or s P_K.
- (L2) Every entry of Pre_K(t) is 0 or a single signed monomial.

Base case: a leaf gives y += s f x.

Step: Pre_K(t) = Pre_child(t) ∘ T, where T is the composition of the completed phases.

Sign +:
- ψ=U: T = I.
- ψ=V: T is z += Bx. Then Pre[w, x] = Pre_V[w, z] B.
- ψ=U': T is z += Bx followed by y += A z, so row y = (1, A, AB). Then
  Pre[z, x] = Pre_U'[z, x] + Pre_U'[z, z] B, which lies in {-B, 0} + B = {0, B} by (S) for U'.
  Rows below z have Pre_U'[w, z] = 0 by (F1).
- ψ=V': T leaves z = z and sets y = y + A z + AB x. Then Pre[y, z] = A + Pre_V'[y, z], which
  lies in {A, 0} by (S). Also Pre[y, x] = AB.

Sign −:
- ψ=V: T = I.
- ψ=U: row y = (1, A).
- ψ=V': T sets y = y + A z and z = z + Bx. Then Pre[y, x] = Pre_V'[y, z] B, which lies in
  {0, -AB}. Also Pre[y, z] = A + Pre_V'[y, z], which lies in {A, 0}.
- ψ=U': row y = (1, -AB). Then Pre[z, x] = Pre_U'[z, x] + B, which lies in {B, 0}.

In each case (S) holds with P_K = AB.

**Lemma L1 (exit).** Inject at the destination of a leaf inside K. At K's end, every entry has
mass at most 1, and the support lies in {y_K} ∪ scratch(K). The proof is by induction. The
child's exit column has x_K = 0. The remaining phases at K's level then add at most one term to
an entry that was 0, or they cancel:
- Sign +, from U: V adds A c_z to y, and V' then removes it.
- Sign +, from U': y = -A c_z.
- Sign −, from U: y = -A c_z.
- All other cases are unchanged.

**Main bound.** Let K (size m >= 2) be the lowest call that contains both cuts. Let a lie in
phase φ_α and b in a later phase φ_β. Then

M(a, b) = Pre_φβ(b) ∘ T_mid ∘ Suf_φα(a).

By L1 and (F2), the incoming column q at φ_β has every entry of mass at most 1. Its source and
destination entries are as follows.

| K sign | (α, β) | src, dst entries |
|---|---|---|
| + | (U, V), (U', V') | (c_z, 0) |
| + | (U, U') | (0, c_z) |
| + | (U, V') | (c_z, A c_z) |
| + | (V, U') | (0, 0) |
| + | (V, V') | (0, c_y) |
| − | (V, U), (V, U'), (V', U') | (0, 0) |
| − | (V, V') | (0, c_y) |
| − | (U, V') | (c_z, 0) |
| − | (U, U') | (0, c_z) |

Registers outside φ_β's footprint keep mass at most 1.

By L2, every row satisfies

|row| <= Σ_v |Pre[w, v]| |q_v| <= #(dirty registers in the footprint) <= d_β + [src dirty] + [dst dirty],

where d_β = ceil(log2 |φ_β|) is the number of scratch levels of φ_β.

The related case (U, V') needs one more step. Row y's c_z-part is (A + Pre[y, z]) c_z, which lies
in {A c_z, 0} by (S) for the sign − call V'. Every other row has Pre[w, y] = 0. So every case
has mass at most d_β + 1 <= ceil(log2 ceil(m/2)) + 1 = ceil(log2 m) <= ceil(log2 n). ∎

Not covered: internal leaf intermediates, multiple or accumulated errors, arbitrary joint
errors, rounding, aliasing, the precision invoice, and kappa.

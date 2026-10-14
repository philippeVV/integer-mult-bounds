# expH working notes (progress log; restart-safe)
- d1587.txt saved (ctoth FINDING). expG/MODEL.md not yet present at 01:12.
- Model adopted (to be checked vs expG):
  * ports t in P112 (E8 lines minus frame F0), labels u_t in F2^9, u_s.u_t=(<r_s,r_t>+1) mod 2.
  * scalar req: net centre matrix M (=D T^T up to scale) has M[S,S]=1, M[S,T]=0 for every pair NOT side-wire-able.
    side wire T->S legal iff source label subset target cap: line: u_T in u_S^perp (u.u=0 <=> odd roots).
    plane: E_T subset E_S^perp (full plane orthogonality). So free pairs shrink.
  * term cost = rank of form on U_k x V_k, U_k=sum E_T (row supp), V_k=sum E_S (col supp)  [#1516 floor; Thm C equality in line case].
  * diagonal-covering term costs >= rank(E_T x E_T)=2 (E_T nondeg). naive private = 112*2 = 224.
- Frame decoder: centre k <-> frame line f_k; row/col supp = 56 ports with <r,f_k> odd; U_k = u_{f_k}^perp (dim 8, rank 8).
  Idea A: r_T in cap_{k in supp T} u_{f_k}^perp keeps the 8 centres at 64; then correct broken odd pairs.
- 01:25 e8.py: labels via E8/2E8 (u=(x,1), form b+1), general Gram G (not std dot). Baseline frame decoder reproduced: 8x56 supports, cost 64, 2G exact.
- Unified model: side wire = single-entry term, cost rank(E_T x E_S) (0 iff planes orth). Decoder = terms summing to EXACTLY I.
- Constant r (rho = perp of frame labels): every term costs line-cost+1, no free pairs -> >=112 terms -> C>=176-ish. FAIL.
- KEY IDEA (tensor 2-colouring): d_t = r_t (x) a_t in C^8 (x) C^2, a_t in {e1,e2} (class). 16 terms, each cost<=9 => C<=144<163.6
  PROVIDED all odd pairs inside each class are plane-orthogonal (E_S perp E_T). Cross-class odd pairs are cancelled by tensor (a_S.a_T=0).
  -> CSP: find class bit + r_t per port. Searching.
- 01:40 SA (sa.py) can't reach 0; per port, best r keeps <=28 of 52 odd nbrs plane-orth (24/26/28 options).
- sat2.py (pycryptosat via uv py3.12): 2-class tensor exact on all 112 = UNSAT (fast). [encoding sanity check still TODO]
- Next: drop set P (private diag terms cost 2 each): C = 2|P| + sum(16 tensor term costs). Need |P|<=9 (cost-9 terms) or <=17 (cost-8). sat3.py with cardinality.
- 02:10 sat3 m=8 UNSAT. m=16/24/32 running (bg task bw22ckbp6, slow). sa_drop weak (drop~78).
- verify_tensor.py independent checker works (bad=0 on SA sols). tensor term costs 5-9.
- price.py (first-order, non-centre word fixed): C=144(16x9) -> a_cap x1.54; 128 -> x1.98; 160 -> x1.09; 162 -> x1.036; 224 rejected.
- 03:00 crash (OOM) killed sat3 runs; now use ../capped, <=4 workers. odd graph: 25088 K3, 88704 K4, 161280 K5 (k5.json). k5col.py = necessary 2-colour cond.
- maxfree.py: best free odd edges 1120 (unrestricted r), 896 (r in W_t). tensor in-class odd edges ~1400+ needed.
- 03:40 k5col: K5-free 2-colouring of odd graph UNSAT at m=0 (clean proof tensor-2class impossible: >=5 mutually orth planes impossible in F2^9). m>=8 timed out.
- odd graph eigen 52,8,4,-4 -> maxcut<=1568 -> in-class odd >=1344 edges for any 2-partition.
- sa2.py (good SA): tensor+drops needs ~56 drops -> C ~ 128+112 = 240 > 224. Tensor route DEAD.
- frame(64)+row-star corrections >= 64+112 = 176 > 163.6 (each port needs >=1 cost-1 correction unless multi-row terms).
- 04:30 frame+row-stars (framestar.py, exact, bad=0): C=288 (r in W_t: frame 64 + 112 stars of cost 2); unrestricted 294.
  König: stars needed >=112 => frame+stars >= 176 > 163.6 regardless.
- starsearch.py: could not get any star to cost 1 under r in W_t (best 224 stars).
- single class (sa1c.py): |Q|<=32, C=228..232. two-class tensor+drops (sa2.py): 240..264. All verified exact by verify_tensor.py.
- minrank float probe (minrank_probe.py, best-free r=mf_1_0, 2240 free entries): rank 90/100 completions exist (res 1e-13), K=70 res 0.086, K=40 res 4.1. => ~70-90 terms (float only).
- BEST C FOUND = 224 (naive). Nothing below.

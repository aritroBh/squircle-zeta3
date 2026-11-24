work = "<project>/research/08_theory/TOW1R_review/"
report = r'''# TOW1R — adversarial referee report: the all-k proof of the weight-three law for M_k

## OVERALL VERDICT: PROOF COMPLETE.

I found **no mathematical error**. The chain
  tower law (TOW1, depths δ ≤ 2, all d, all k in the stated ranges)
    ⟹ affine Gm law for the degree-3 keys (3,0),(2,1), with slopes s30,s21
    ⟹ WT3's reduction theorem (identities C1–C4)
    ⟹ the kernel law R_20(k;v), R_11(k;v) and W3(k) = W3' − (k−3)(ζ(3)/8 + π²log2/24 + log³2/12) for all k ≥ 3
is logically sound, and every "PROVED for all d" step genuinely is (not a hidden finite check).
I independently confirmed the *conclusion* on data OUTSIDE every range the authors used:
exact degree-2 kernels at **k = 16, 18, 20**; the affine value law by 58-digit quadrature at
**k = 15, 16, 18, 20**; the tower law and the affine Gm law at **k = 15, 16, 17, 18, 20**.
The only issues are wording/documentation (stale status lines, one wrong identity in a remark). See §W.

## Item 1 — logic chain & slope determination: OK
- k-ranges are consistent. The theorem needs k ≥ 3; both degree-3 keys (3,0) and (2,1) are
  off-diagonal, so their tower law runs from k = d = 3. (2,1) is adjacent; its diagonal parent
  (2,2) is exceptional only at k = 4 (=2p), handled by blindness E3 (p=2). No k ≥ 3 is missed.
- **The slopes s30, s21 are NOT fit to finite data.** The tower law gives u(δ=0,1) constant and
  u(δ=2) affine in k, with e(2,1) = −ρ(0)/(2d−3) (proved: the (δ,s)=(2,1) Newton equation E2a
  reduces to E0, WT2b's Lemma D). Hence Gm_30,Gm_21 are affine in k, and the slope is a single
  k-independent function determined by the closed forms. I DERIVED it myself:
  slope of Gm_30 = (√2/12)(1/(v−i)+1/(v+i)) − 1/12 = −(1/12)(1 − 2√2 v/(v²+1)) = s30, and
  s21 = −3 s30 (verified as PFs, step_15; and the affine law Gm30^(k)=Gm30^(3)+(k−3)s30 checked
  EXACTLY at k=18,20, step_26). The write-up says "verified exactly k≤14"; this is only a
  confirmation — the value is derivable from proved quantities (recommend stating the derivation).
- The reduction theorem P3: R_20(k)=R_20(3)+(k−3)S_20, R_11(k)=R_11(3)+(k−3)S_11 follows from the
  affine Gm law + C1–C4 (all verified independently, step_15). The value slope integral
  ∫(S_20 log²t + S_11 log t log(1−t))dv = −ζ(3)/8 − π²log2/24 − log³2/12 confirmed to 58 digits (step_03).
- **No circularity**: E0 (WT2b Lemma D/I1–I4), Lemma C (ρ(0)), Theorem ρ (ρ(1),ρ(2) via
  log-derivatives), E2a(=E0), E1, E2b, E3 are all proved without assuming the tower law or WT3.

## Item 2 — the symbolic all-d proofs (E1, E2b, E3): OK
- The scripts treat d,p as genuine sympy symbols; a remainder that simplifies to 0 is a complete proof.
- I cross-checked that each script encodes EXACTLY the Newton identity the proof needs: I rebuilt
  every script's LHS and RHS, evaluated at concrete (d,p), and matched my OWN independent engine's
  e-coefficients / Λ-feeds (normalized) — E1-adj (step_11), E2b-nonadj & E2b-adj (step_12), E3 (step_14).
- I re-derived E1 non-adjacent by hand AND symbolically (→0, step_07). I re-verified ALL the Newton
  identities E0,E1,E2a,E2b with my own encoding for **d = 2..25** (authors: d ≤ 11) — 669 cells, 0 bad
  (step_09). E3 blindness feed = 0 with my own encoding for p = 2..12 (step_13). Building blocks
  (ρ-ratios rA0/rB0/E3-ratio, τ-values, conj relations) verified exactly (steps 06, 10).

## Item 3 — initial data ρ(0),ρ(1),ρ(2) and S(t*),S'(t*): OK
- ρ(0),ρ(1),ρ(2) match the actual integrand +i principal part (build_integrand_pf) for d = 2..12,
  137 cells, 0 bad (step_05). ρ(1) = −i(d−2)ρ(0) and the ρ(2) closed form confirmed.
- Corollary S3: (S′/S)(t*) = r_m = 2(d−2m) − 2di − 4(d−1)(d−2m)/(2d−3) verified d = 1..14 (step_18).
- The log-derivative derivation is sound: t′(i)=0, t″(i)=−i/2, 1/(1+a²)+1/(1+b²)=1, the
  −d/dv[…](i) = −1 term, and A″/A = (3d−3)/2 − (i/2)r_m − (d−2)² all check out (steps 18,20,19-pieces).
- WT2b Lemma C (ρ(0)) and the stabilization identities I1,I2,I3,I4 (= E0) verified for p,d ≤ 29
  (steps 20,22), including the CX4 intermediate relations and ρ_pp pure-imaginary (step_27).

## Item 4 — fold identity & E3 (edge cases): OK
- Fold lemma DL1(1/v)/v² = −DL2(v) PROVED exactly (step_07); hence the §6bis identity
  −p·DL2·[Gm(v)−Gm(1/v)]. The T1 fold formula Fold_p matches the actual reducer fold (pf_fold)
  for k = 4..12, 107 cells, 0 bad (step_16).
- E3 (diagonal blindness) proved all p (script), cross-checked against my numeric feed (step_14),
  and holds with my own encoding p = 2..12 (step_13). p = 1 is correctly out of scope (the child
  (1,0) has degree 1, no HigherPart; E3 concerns adjacent keys (p,p−1), p ≥ 2). The binom(4p−6,2)→
  (4p−6)(4p−7)/2 replacement is valid for p ≥ 2.
- First diagonal levels verified from data: diagonal base case u_pp^(2p+1)=ê(δ,0); diagonal tower
  law k = 2p+1..2p+4; blindness of the adjacent child at k = 2p — all 0 bad (step_21).

## Item 5 — T1 uses CX4's corrected Lemma B: OK
- At depth 0, T1's direct-DL2 term (indicator [p>q]) plus the fold combine to p·τ2·(u_pp − conj u_pp),
  i.e. exactly CX4's corrected numerator with the [i ≥ j+2] indicator (verified by hand).
- The full T1 recursion (DL1 + DL2-direct + fold), my OWN implementation, matches tower data for
  k = 3..12, 1580 cells, 0 bad (step_17). No double-counting at adjacent keys.

## Item 6 — independent data check (outside all author ranges): OK
- Degree-2 kernels R_20(k;v), R_11(k;v) match the claimed closed forms EXACTLY at **k = 16, 18, 20**
  (exact Q(√2,i) arithmetic; higher part empty) — step_02. Authors checked only k ≤ 15.
- Value form: W3(k) − W3(3) = (k−3)·slope to ~60 digits at k = 15,16,18,20; slope integral
  = −ζ(3)/8 − π²log2/24 − log³2/12 to 58 digits — step_03.
- Tower law U_e == data at k = 15,16,17 (621 cells, step_08) and degree-3 tower law at k = 18,20
  (step_25); affine Gm law Gm30^(k)=Gm30^(3)+(k−3)s30 EXACT at k = 18,20 (step_26).
- Control: the fast reducer equals the paper's own ibp_reduce.reduce_k on ALL keys at k = 3,4
  (step_23), so the machinery is faithful. Authors' verify_tower.py 12/12 and verify_wt3.py 10/10
  both rerun and PASS (steps 04, 24).

## W — wording / documentation fixes (none affects the math)
1. proof_tower.md §1: "τ2_n = conj(τ1_n)" is WRONG for odd n. Correct: conj(τ1_n) = (−1)ⁿ τ2_n
   (equivalently τ2_n = (−1)ⁿ conj(τ1_n); follows from DL2(v) = −DL1(−v)). Verified numerically.
   The code uses the correct relation (tc1c); only the prose is wrong.
2. proof_tower.md status header: "E1, E2b, E3 are currently VERIFIED for d ≤ 11" is STALE —
   §6 and prove_*.py PROVE them for all d. Update the header.
3. proof_tower.md §7 / result.md: the slopes are described as "verified exactly"; they are in fact
   DERIVABLE in closed form from the tower law (e(2,1) = −ρ(0)/(2d−3)) and ρ(0). State the derivation
   (one line: slope = −e(2,1)-feed = ρ(0)/(2d−3) integrated and normalized).
4. WT3 result.md / proof_wt3.md: still call the affine Gm law "the single remaining lemma /
   REMAINING". Now that TOW1 proves it, update WT3's status to "PROVED (via TOW1)".
5. proof_tower.md "reduction to identities" theorem: the hypothesis list omits (E2a) though the
   proof uses it ("the slope equation at s = 1"); it is proved free (= E0) immediately after.
   Harmless, but list (E2a) explicitly for clarity.

## Evidence table (all exact arithmetic unless a digit count is stated; residuals 0 unless stated)
| check | range | result | my step |
|---|---|---|---|
| kernels R_20,R_11 vs claimed | k=3,4,5 | exact match | step_01 |
| **kernels R_20,R_11 vs claimed (NEW DATA)** | k=16,18,20 | exact match | step_02 |
| W3(k) affine + slope integral (quadrature) | k=15,16,18,20 | res ≤ 6e-61 (60 dps) | step_03 |
| authors' verify_tower.py rerun | — | 12/12 PASS | step_04 |
| ρ(0),ρ(1),ρ(2) vs integrand | d=2..12 | 137 cells, 0 bad | step_05 |
| τ1_n,τ2_n, t-dict, conj(τ1_n)=(−1)ⁿτ2_n | n≤3 | match | step_06 |
| fold lemma DL1(1/v)/v²=−DL2; E1 non-adj | exact | 0 / proved | step_07 |
| **tower law U_e==data (NEW k)** | k=15,16,17 | 621 cells, 0 bad | step_08 |
| Newton identities E0,E1,E2a,E2b (my encoding) | d=2..25 | 669 cells, 0 bad | step_09 |
| ρ-ratios rA0,rB0,E3-ratio, conj purity | p,d≤29 | 0 bad | step_10 |
| prove_E1_adj encodes correct identity | p=3,4,5 | match | step_11 |
| prove_E2b (non)adj encode correct identity | d,p | match | step_12 |
| E3 blindness feed=0 (my encoding) | p=2..12 | 0 bad | step_13 |
| prove_E3 encodes correct feed | p=3,4,5 | match (=0) | step_14 |
| C1,C2,C3,C4 + s30 PF derivation | exact | all True | step_15 |
| fold formula vs reducer fold | k=4..12 | 107 cells, 0 bad | step_16 |
| FULL T1 recursion (my impl) | k=3..12 | 1580 cells, 0 bad | step_17 |
| Corollary S3 (S′/S)(t*)=r_m | d=1..14 | 0 bad | step_18 |
| ρ(1),ρ(2) derivation pieces (t′,t″,1/(1+a²)+…,−1,I1,sum) | exact | verified | step_20 |
| diagonal base/tower-law/blindness (incl p=1) | p≤5 | 0 bad | step_21 |
| WT2b I2,I3,I4 stabilization | p,d≤29 | 0 bad | step_22 |
| CONTROL fastred == paper ibp_reduce, ALL keys | k=3,4 | 0 bad | step_23 |
| authors' verify_wt3.py rerun | — | 10/10 PASS | step_24 |
| degree-3 tower law δ≤2 (NEW k) | k=18,20 | 0 bad | step_25 |
| **affine Gm law Gm30 (NEW k)** | k=18,20 | EXACT | step_26 |
| CX4 I3/I4 intermediate relations | p≤24 | 0 bad | step_27 |

## Scripts
myengine.py (my independent e/ê/Λ/ρ/τ engine); step_01..step_27.py (each self-contained; the
harness also auto-saved my inline blocks as step_28+). All ran under nice, one at a time, ≤ 30 min.
'''
with open(work+"result.md","w") as f: f.write(report)
print("result.md written,", len(report), "chars")
with open(work+"progress.log","a") as fh: fh.write("step_28: wrote result.md (verdict PROOF COMPLETE)\n")

# TOW1R — adversarial referee report: the all‑k proof of the weight‑three law for M_k

## OVERALL VERDICT: PROOF COMPLETE.

I found **no mathematical error**. The chain
> tower law (TOW1, depths δ ≤ 2, all d, all k in the stated ranges) ⟹ affine Gm law for the degree‑3 keys (3,0),(2,1) with slopes s30,s21 ⟹ WT3's reduction theorem (identities C1–C4) ⟹ the kernel law R_20(k;v), R_11(k;v) and W3(k) = W3' − (k−3)(ζ(3)/8 + π²log2/24 + log³2/12) for all k ≥ 3

is logically sound, and every "PROVED for all d" step genuinely is (not a hidden finite check). I independently confirmed the **conclusion on data outside every range the authors used**: exact degree‑2 kernels at **k = 16, 18, 20** (authors: k ≤ 15); the affine value law by 60–80‑digit quadrature at **k = 15, 16, 17, 18, 19, 20**; the tower law and the affine Gm law at **k = 15, 16, 17, 18, 20** (authors: k ≤ 14). The only issues are wording/documentation (§W); none affects the mathematics.

## Item 1 — logic chain & slope determination: OK
- **k‑ranges consistent.** Theorem needs k ≥ 3; both degree‑3 keys are off‑diagonal so their tower law runs from k = d = 3. (2,1) is adjacent; its diagonal parent (2,2) is exceptional only at k = 4 (=2p), handled by blindness E3 (p=2). No k ≥ 3 is missed.
- **Slopes s30, s21 are NOT fit to finite data** (the point flagged in the brief). The tower law makes u(δ=0,1) constant and u(δ=2) affine in k, with e(2,1) = −ρ(0)/(2d−3) (proved: the (δ,s)=(2,1) Newton equation E2a reduces to E0 = WT2b Lemma D). Hence Gm_30, Gm_21 are affine with a single k‑independent slope function determined by closed forms. I **derived** it: slope(Gm_30) = (√2/12)(1/(v−i)+1/(v+i)) − 1/12 = −(1/12)(1 − 2√2 v/(v²+1)) = s30, and s21 = −3s30 (PF check, step_15; affine law Gm30^(k)=Gm30^(3)+(k−3)s30 EXACT at k=18,20, step_26). The write‑up's "verified exactly k≤14" is only a confirmation — the value is derivable from proved quantities.
- Reduction theorem P3 (R_ij(k) = R_ij(3) + (k−3)S_ij) follows from the affine law + C1–C4 (all verified, step_15). Value slope integral ∫(S_20 log²t + S_11 log t log(1−t))dv = −ζ(3)/8 − π²log2/24 − log³2/12 confirmed to 77 digits (step_29).
- **No circularity**: E0 (WT2b I1–I4), Lemma C (ρ(0)), Theorem ρ (ρ(1),ρ(2) via log‑derivatives), E2a(=E0), E1, E2b, E3 are all proved without assuming the tower law or WT3.

## Item 2 — symbolic all‑d proofs (E1, E2b, E3): OK
Scripts treat d,p as genuine symbols; remainder→0 is a complete proof. I rebuilt each script's LHS/RHS, evaluated at concrete (d,p), and matched my **own independent engine's** e/ê/Λ‑feeds (normalized): E1‑adj (step_11), E2b‑nonadj & E2b‑adj (step_12), E3 (step_14) — all match, so the scripts prove the *right* identities. I re‑derived E1 non‑adjacent by hand and symbolically (→0, step_07); re‑verified **all** Newton identities E0,E1,E2a,E2b with my own encoding for **d = 2..25** (authors: d ≤ 11), 669 cells 0 bad (step_09); E3 feed=0 with my encoding p = 2..12 (step_13). Building blocks (ρ‑ratios, τ‑values, conj relations) verified exactly (steps 06, 10).

## Item 3 — initial data ρ(0),ρ(1),ρ(2), S(t*), S'(t*): OK
ρ(0),ρ(1),ρ(2) match the actual integrand +i principal part (build_integrand_pf) for d = 2..12, 137 cells 0 bad (step_05). Corollary S3 (S′/S)(t*) = r_m verified d = 1..14 (step_18). Log‑derivative derivation pieces all check: t′(i)=0, t″(i)=−i/2, 1/(1+a²)+1/(1+b²)=1, the −1 term, A″/A = (3d−3)/2 − (i/2)r_m − (d−2)² (steps 18,20). WT2b Lemma C (ρ(0)) and stabilization I1–I4 (=E0) verified for p,d ≤ 29, incl. CX4's intermediate relations and ρ_pp pure‑imaginary (steps 20,22,27).

## Item 4 — fold identity & E3 (edge cases): OK
Fold lemma DL1(1/v)/v² = −DL2(v) PROVED (step_07) ⟹ §6bis identity −p·DL2·[Gm(v)−Gm(1/v)]. Fold formula Fold_p matches the actual reducer fold for k = 4..12, 107 cells 0 bad (step_16). E3 proved all p (script), cross‑checked vs my numeric feed (step_14), holds with my encoding p = 2..12 (step_13). p = 1 is correctly out of scope (child (1,0) has degree 1, no HigherPart); the binom(4p−6,2)→(4p−6)(4p−7)/2 replacement is valid for p ≥ 2. First diagonal levels verified from data: base case u_pp^(2p+1)=ê(δ,0), diagonal tower law k = 2p+1..2p+4, blindness of adjacent child at k = 2p — all 0 bad (step_21).

## Item 5 — T1 uses CX4's corrected Lemma B: OK
At depth 0, T1's direct‑DL2 term (indicator [p>q]) plus the fold combine to p·τ2·(u_pp − conj u_pp), exactly CX4's corrected numerator with the [i ≥ j+2] indicator (verified by hand). Full T1 recursion (DL1 + DL2‑direct + fold), my **own** implementation, matches tower data for k = 3..12, 1580 cells 0 bad (step_17). No double‑counting at adjacent keys.

## Item 6 — independent data check (outside all author ranges): OK
Degree‑2 kernels R_20(k;v), R_11(k;v) match the claimed closed forms **EXACTLY at k = 16, 18, 20** (exact Q(√2,i); higher part empty) — step_02. Value form W3(k)−W3(3) = (k−3)·slope to ~80 digits and slope integral to 77 digits at k = 15..20 — steps 03, 29. Tower law U_e==data at k = 15,16,17 (621 cells, step_08) and degree‑3 tower law at k = 18,20 (step_25); affine Gm law EXACT at k = 18,20 (step_26). Control: fast reducer == paper's own ibp_reduce.reduce_k on ALL keys at k = 3,4 (step_23), so the machinery is faithful. Authors' verify_tower.py 12/12 and verify_wt3.py 10/10 rerun and PASS (steps 04, 24; the earlier claude_rerun_verify4.out was truncated, and the "two M1 fixes" are an import and a stale‑k variable — both harmless verifier‑script fixes, not proof changes).

## W — wording / documentation fixes (none affects the math)
1. **proof_tower.md §1**: "τ2_n = conj(τ1_n)" is WRONG for odd n. Correct: conj(τ1_n) = (−1)ⁿ τ2_n (equivalently τ2_n = (−1)ⁿ conj(τ1_n); follows from DL2(v) = −DL1(−v)). Verified numerically; the code uses the correct relation (tc1c).
2. **proof_tower.md status header**: "E1, E2b, E3 currently VERIFIED for d ≤ 11" is STALE — §6 and prove_*.py PROVE them for all d. Update the header.
3. **proof_tower.md §7 / result.md**: slopes described as "verified exactly"; they are DERIVABLE in closed form (slope = ρ(0)/(2d−3) integrated and normalized). State the one‑line derivation.
4. **WT3 result.md / proof_wt3.md**: still call the affine Gm law "the single remaining lemma / REMAINING". Now that TOW1 proves it, update to "PROVED (via TOW1)".
5. **proof_tower.md "reduction to identities" theorem**: the hypothesis list omits (E2a) though the proof uses it ("the slope equation at s = 1"); it is proved free (= E0) immediately after. List (E2a) explicitly for clarity.

## Evidence table (exact arithmetic unless a digit count is stated; residuals 0 unless stated)
| check | range | result | step |
|---|---|---|---|
| kernels R_20,R_11 vs claimed | k=3,4,5 | exact match | step_01 |
| **kernels R_20,R_11 vs claimed (NEW DATA)** | k=16,18,20 | exact match | step_02 |
| W3(k) affine + slope integral (quadrature) | k=15,16,18,20 | res ≤ 6e‑61 (60 dps) | step_03 |
| authors' verify_tower.py rerun | — | 12/12 PASS | step_04 |
| ρ(0),ρ(1),ρ(2) vs integrand | d=2..12 | 137 cells, 0 bad | step_05 |
| τ1_n,τ2_n, t‑dict, conj(τ1_n)=(−1)ⁿτ2_n | n≤3 | match | step_06 |
| fold lemma DL1(1/v)/v²=−DL2; E1 non‑adj | exact | proved / 0 | step_07 |
| **tower law U_e==data (NEW k)** | k=15,16,17 | 621 cells, 0 bad | step_08 |
| Newton identities E0,E1,E2a,E2b (my encoding) | d=2..25 | 669 cells, 0 bad | step_09 |
| ρ‑ratios rA0,rB0,E3‑ratio, conj purity | p,d≤29 | 0 bad | step_10 |
| prove_E1_adj encodes correct identity | p=3,4,5 | match | step_11 |
| prove_E2b (non)adj encode correct identity | d,p | match | step_12 |
| E3 blindness feed=0 (my encoding) | p=2..12 | 0 bad | step_13 |
| prove_E3 encodes correct feed | p=3,4,5 | match (=0) | step_14 |
| C1,C2,C3,C4 + s30 PF derivation | exact | all True | step_15 |
| fold formula vs reducer fold | k=4..12 | 107 cells, 0 bad | step_16 |
| FULL T1 recursion (my impl) | k=3..12 | 1580 cells, 0 bad | step_17 |
| Corollary S3 (S′/S)(t*)=r_m | d=1..14 | 0 bad | step_18 |
| ρ(1),ρ(2) derivation pieces (t′,t″,1/(1+a²)+1/(1+b²),−1,I1,sum) | exact | verified | step_20 |
| diagonal base/tower‑law/blindness (incl p=1) | p≤5 | 0 bad | step_21 |
| WT2b I2,I3,I4 stabilization | p,d≤29 | 0 bad | step_22 |
| CONTROL fastred == paper ibp_reduce, ALL keys | k=3,4 | 0 bad | step_23 |
| authors' verify_wt3.py rerun | — | 10/10 PASS | step_24 |
| degree‑3 tower law δ≤2 (NEW k) | k=18,20 | 0 bad | step_25 |
| **affine Gm law Gm30 (NEW k)** | k=18,20 | EXACT | step_26 |
| CX4 I3/I4 intermediate relations | p≤24 | 0 bad | step_27 |
| W3(k) affine + slope integral (quadrature, NEW k) | k=17,19 | res ≤ 9e‑81 (80 dps); slope res 1.4e‑77 | step_29 |

## Scripts (deliverables, all in new_results/TOW1R_review/)
- **myengine.py** — my independent e/ê/Λ/ρ/τ engine (closed forms + Newton recursion + fold), built from first principles and validated against data.
- **step_01..step_27.py, step_29.py** — each self‑contained, run under `nice`, one at a time, ≤ 30 min: step_01/02 (kernels), step_03/29 (quadrature), step_04/24 (author verifier reruns), step_05 (initial data), step_06 (τ), step_07 (fold lemma + E1), step_08 (tower law new k), step_09 (Newton d≤25), step_10 (ρ‑ratios), step_11/12/14 (prove‑script cross‑checks), step_13 (E3), step_15 (C1–C4 + slopes), step_16/17 (fold + full T1), step_18/20 (S3 + ρ derivation), step_21 (diagonal edge cases), step_22/27 (I2/I3/I4 + relations), step_23 (paper control), step_25/26 (new‑k tower + affine Gm).
- The harness also auto‑saved my inline code blocks as step_28+ files.

**Bottom line:** the all‑k proof of the weight‑three law is correct and complete. The degree‑2 kernels R_20(k;v), R_11(k;v) and W3(k) = W3' − (k−3)(ζ(3)/8 + π²log2/24 + log³2/12) are proved for every k ≥ 3. The five wording items in §W should be applied to the write‑ups; they do not affect the validity of the proof.
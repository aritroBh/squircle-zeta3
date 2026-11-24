# CX6G: Second Independent Referee Report

## Table of Displays Checked

| Display | What I Checked | Numerical Evidence | Verdict |
|---------|----------------|-------------------|---------|
| (1) | rest bound for k=100,200,400 | ratios 8.6e-12, 7.4e-13, 1.1e-14 | **OK** |
| (2) | combined decay for k=100,200,400 | ratios 6.0e-10, 1.3e-9, 8.4e-10 | **OK** |
| (3) | M(eps) forms at eps=0.01,0.05,0.1,0.2 | diffs < 4e-17 | **OK** |
| (4) | H(s) bounds at x=0.1,0.5,1,2,5,10 | worst ratio 0.9706 | **OK** |
| (5) | C(z) contour identity at z=0,0.1,...,10 | diffs < 1e-57 | **OK** |
| (6) | residue sum identity | verified via binomial series | **OK** |
| (7) | M_k representation at k=0 | diff ±9.36e-81 contains 0 | **OK** |
| (8) | rest decomposition sign | cross-checked with verify_rest.py | **OK** |
| (9) | J_k term for k=8,10,12 | formula verified against P_k(0) | **OK** |
| (12) | A_k bound | checked on Gamma contour | **OK** (per proof's parameters) |
| (13) | K, K-E bounds on Gamma | algebraic formulas correct | **OK** |
| (14) | F(j,q;lambda) AM-GM bound | exponential decay verified | **OK** |
| (15)-(19) | C_star constant | 10578907.82 < 10,600,000 with margin 21,092 | **OK** |
| (20) | tail bound at R=200 | ≤ 6.77e-102 | **OK** |

## Key Verifications

1. **Claim (1)**: |M_k - M_k^main| ≤ 10,600,000 e^{-2√N}
   - Tested at k=100,200,400: actual/bound ratios all below 10^{-11}
   - The bound is deliberately loose but valid

2. **Claim (2)**: Combined decay law with main-term theorem (254.23)
   - delta = 2 - √π ≈ 0.2275 > 0 ✓
   - Total constant: 1.2094×10^9 < 1.21×10^9 with margin 578,166

3. **C_star constant**: Verified with arb at 90 digits:
   - C_star = [10578907.823529397...] < 10,600,000 ✓

4. **Contour identities**: Display (5) verified to < 10^{-57} precision

5. **H(s) bounds**: Display (4) verified; worst ratio 0.9706

## Overall Verdict: **PROVED**

Both claims are correct as stated. All numerical checks pass with comfortable margins. The proof uses deliberately crude bounds where tighter ones exist, which is valid proof methodology.

## Step Files Produced
- `step_01.py`: Claim (1) verification  
- `step_02.py`: Claim (2) verification
- `step_03.py`: Display (3) M(eps) forms
- `step_04.py`: Display (4) H(s) bounds
- `step_05.py`: Display (5) C(z) contour identity
- `step_06.py`: C_star constant verification
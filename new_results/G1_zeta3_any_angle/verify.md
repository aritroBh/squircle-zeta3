# G1 verification: ζ(3) at a corner of any angle

**Verdict: CONFIRMED WITH CORRECTIONS.**
- Every number in the solver's result.md that I re-tested holds.
- The PROVED headline claims hold. Theorem A, the canonical ζ(3) coefficient (0 in M₂, −Φ(c) in M₃), and Theorem C each
  survived a fresh derivation or a check by a different method.
- The corrections are about rigour and wording. None changes a number.
- Two items can be **upgraded**: B1 is now proved for all n (§3.4), and W₂ has a 5-line proof (§3.3).

Working folder: `<scratch>/verify_G1_zeta3_any_angle/`. It holds a copy of the solver's folder plus my
`v_*.py` scripts and their `.out` files. Durable copies are in `G1_zeta3_any_angle/verify/`. Every number below was
printed by one of those scripts, and the file is named in each line.

---

## 1. Re-runs of the solver's scripts (copied to the scratch folder)

| Script | Printed result | Matches result.md? |
|---|---|---|
| `m2_structure.py` (`rerun_m2_structure.out`) | key identity, IBP identity (V·(log g)′ − 1/(2(t+1)) = 0, dV/dt − R_A^c = 0), g(t₁) = 1; residues of R_A^s, R_B, R_C as tabulated; partial-fraction remainders 0; P does not divide den(R_A^s) | yes |
| `m2_w2_identity.py` | log part S(M+P) = log s·log c²; numeric ≤ 1.9e-121 at 9 angles | yes |
| `flat_symbolic_proof.py 3` | ζ(2)-coeff of μ₂,ₙ = 1/((2n−1)(2n+1)); ζ(3)-coeff of μ₃,ₙ = −1/((2n−1)(2n+1)) | yes |
| `flat_mu2_symbolic.py` | harmonic-number formula = exact table for all n ≤ 320: True | yes (see note) |
| `flat_sum_check.py` | partial sums vs series: 4.5e-322 (5π/6) … 2.3e-102 (π/2) | yes |
| `m1_classical_check.py` | ≤ 5.1e-321 at 9 angles | yes |
| `zeta3_bookkeeping.py` | Z(c) = (8c+7)/4; flat-normalised [0, 0, 0]; squircle −19/16 + 189√2/64; 2π/3: 3/16 | yes |
| `report_formula_check.py` | B1, B2′, 2π/3: ≤ 1.1e-120; HPL B2: ≤ 1.0e-61 | yes |
| `m2_2pi3_pslq.py` (300 digits) | bases A/B/C: max coefficient 408/369/180, residual 2.15e-321/2.15e-321/1.07e-321; ζ(3) coefficient 3/16, 123/16, −29/4 | yes |

Note: the solver's archived `flat_mu2_symbolic.out` also contains stale "MISMATCH … False" lines from an earlier run.
A clean re-run prints only `True`.

## 2. Independent checks (different method or parametrisation, 46 or more digits)

**2.1 Values by a new method** (`v_cauchy.py`, `v_compare_values.py`).
- Method: a Cauchy contour integral in u of the full corner integral M(u;α), using the form
  (1+h)√(1 − 4c²h/(1+h)²) so that the square root stays analytic. It does not use the solver's u-series-at-nodes code.
- Two runs: (R, N, dps) = (0.4, 128, 55) and (0.33, 160, 60). The runs agree to ≤ 2.3e-40.
- Run B agrees with the solver's 318-digit values for M₀…M₃ at all 9 angles to **≤ 2.0e-46**.
- M₂(π/2) agrees with the squircle 1200-digit value to 4.4e-47.
- Four **fresh angles** (0.3, 0.9, 2.6 and 3.0 rad) were added as test points. c = 0.98877 at 0.3 rad.
- Solver values vs D10b's 340-digit file (`v_d10b_compare.py`): |ΔM₂| = 1.4e-323 (3π/4) and 2.0e-322 (π/6).

**2.2 Theorem A by a different route** (`v_flat_numeric.py`).
- μ_{k,n} was computed numerically from the *original* g-integral on [0,1]. This does not use the y-substitution, the
  g ↦ 1/g symmetry, or Beta derivatives.
- Agreement with the solver's exact table for k ≤ 3, n ≤ 8: **≤ 6.2e-72**. This also confirms the symmetry step and
  the y-form of J_n, which I re-derived by hand.
- PSLQ for n ≤ 3 gives small clean relations with no log 2. Example: μ₃,₁: [−72, 35, −2, −24, 0, 0, 0, 0], i.e.
  35/72 − ζ(2)/36 − ζ(3)/3.
- For n ≥ 4 the 70-digit PSLQ returns spurious relations with coefficients around 10⁷, because the true denominators
  are too large for that precision. They are ignored. The 1e-72 agreement with the exact table is the real check.

**2.3 Theorem A symbolically, own code** (`v_flat_symbolic.py`). I expanded the integrand with explicit binomial series
(the solver used `series(log)`) and reduced it with my own cumulant bookkeeping:
- ζ(2)-coefficient of μ₂,ₙ = **1/((2n−1)(2n+1))**, symbolic in n. All [u²] monomials have log-degree 2.
- ζ(3)-coefficient of μ₃,ₙ = **−1/((2n−1)(2n+1))**, symbolic in n.
- ζ(2)-coefficient of μ₃,ₙ, exact for n ≤ 40, equals both the solver's table and the Taylor coefficients of the solver's
  Ψ(c). The Ψ series has no c^{≤0} or odd terms.
- By hand: the ab-coefficient of [u²] is −n(n−1) + 2n(2n−1)y − 2n(2n−1)y². With B(n+1,n) = B/2 and
  B(n+2,n) = B(n+1)/(2(2n+1)), this gives ζ(2)/(4n²−1).
- Φ(c) = ½ − ((1−c²)/(2c))artanh c, since 1/(4n²−1) = ½[1/(2n−1) − 1/(2n+1)]. Numerically ≤ 4.9e-62 at 13 angles.

**2.4 Closed forms typed by hand from result.md** (`v_closed_forms.py`). Compared with my Cauchy values at 13 angles
(9 plus 4 fresh); worst residuals:

| Formula | worst residual (13 angles) |
|---|---|
| B0 | 4.2e-48 |
| B1 | 1.6e-47 |
| B2′ (uniform classical, Z = 7/4 + 2c) | 5.7e-47 |
| B2 (HPL form, my own HPL power series; S(−c) evaluated directly) | 5.7e-47 |
| W₂ + W₃ (W₃ by quadrature of the x-table of Thm C) | 5.7e-47 |
| W₂ closed form vs its t-integral | 1.2e-61 |

These residuals sit at the noise floor of my Cauchy values. At 330 digits against the solver's 318-digit values:
- the 2π/3 formula agrees to 4.5e-323;
- B2′ at c = √3/2 (π/3) agrees to 4.8e-322;
- B2′ at the other 7 angles agrees to ≤ 3.4e-322, and B1 to ≤ 4.2e-321.

At the fresh angle 0.3 rad (c = 0.98877) we have c⁴⁰⁰ = 0.011 and c⁶⁴⁰ = 7.3e-4 (`v_tail_sizes.out`). So this angle
tests Taylor coefficients far beyond the fitted range, and B2 still agrees to 5e-47.

**2.5 Exact Taylor check of B1/B2, own code** (`v_taylor_exact.py 300`).
- In exact rationals, the Taylor coefficients of the HPL forms of B1 and of the rational part of B2 equal the
  Theorem A harmonic-number formulas (typed from result.md) **for every n ≤ 300**. The solver claimed n ≤ 200.
- There are no odd or negative powers. M₀ also matches −1/(n(2n−1)).
- A mutation test (changing one coefficient by 1/6) makes the check fail, as it should.

**2.6 M₃** (`v_m3_newangles.py`, the solver's 74-term R₃ and Ψ). Against my independent M₃ at 5 angles not used by
the solver (2.0, 0.3, 0.9, 2.6, 3.0 rad): **≤ 2.0e-46**. This includes c = 0.98877, where c⁶⁴⁰ = 7.3e-4. Separately,
Σ_{n≤320} μ_{k,n}c^{2n} from the exact table matches my M₀…M₃ at 7 angles to the noise floor (`v_flat_sum.out`).

**2.7 Basis dependence at 2π/3** (`v_basis_change.py`). The identities that move between bases each hold to 1e-331:
- Landen at x = ⅓;
- inversion at −2;
- duplication at ½;
- Li₃(½).

All three printed forms (ζ(3)-coefficients 3/16, 123/16, −29/4) reproduce M₂(2π/3) to 4.5e-323. By hand:
- trading Li₃(⅓) for Li₃(⅔) adds 15/2, giving 3/16 + 120/16 = 123/16;
- trading Li₃(−½) for Li₃(¼) adds (17/2)(−7/8), giving −29/4.

## 3. Reading the PROVED arguments

**3.1 Theorem A: sound, with one small rigour gap.** I re-derived the y-form of J_n, the identity M(u;π) ≡ 0
(∫₀¹g^{u−1}(1+g)^{−1−u}dg = 2^{−u}/u), and the Beta/polygamma reduction.

*Gap (minor).* Step 2 gives uniform convergence only for real u, using 4c²h/(1+h)² ≤ c². Taking u-Taylor coefficients
term by term needs more: either uniform convergence on a complex disc |u| ≤ δ, or dominated convergence of the
u-derivatives. Either is easy:
- for |u| ≤ δ, |4h/(1+h)²| ≤ 1/(1 − Cδ²), because ρ log²ρ/(1+ρ)² is bounded with ρ = |h|; so c²/(1 − Cδ²) < 1;
- or ∂_u^k brings in only n^k log^k g, and Σ n^k c^{2n} < ∞.

One sentence should be added.

**3.2 The canonical ζ(3) coefficients: correct, but two statements need care.**
- **(Correction 1)** result.md §2 says the M₃ split is unique "given … the standard independence of {1, π², ζ(3)}".
  The ℚ-linear independence of 1, ζ(2), ζ(3) is an **open problem**, not a standard fact. For example, arXiv 1401.5322
  says it "remains an open problem". Rewording needed:
  - The polygamma reduction defines a canonical *formal* expression r′ₙ + q′ₙζ(2) − ζ(3)/(4n²−1). That statement is
    PROVED.
  - Uniqueness of this split for the real number μ₃,ₙ is CONDITIONAL on that conjecture.
  - For M₂ no condition is needed: π² is irrational.
- **(Correction 2: framing)** "Canonical ζ(3)-coefficient of M₂ = 0 for every α" is a statement about the c²-Taylor
  coefficients, i.e. M₂ = ζ(2)Φ + R₂ with R₂ ∈ ℚ[[c²]]. It does not say that the *number* M₂(α) is ζ(3)-free at a
  fixed α. For instance, in the note's basis at π/2 the coefficient is 189√2/64 − 19/16 ≠ 0. result.md explains this in
  §0(c)1; the one-line summary should keep the qualifier "as a power series in c² (flat-angle normalisation)".
- The uniform-basis value Z(c) = 7/4 + 2c is correctly labelled NUMERICAL. Its uniqueness also needs the (open)
  independence of ζ(3), π²log2 and log³2. Its "forcing" argument is valid: Theorem A says no Taylor coefficient
  contains ζ(3), so Z(c) = −Σ C_i(c)z(a_i(0)) as Laurent polynomials. I checked by hand that this gives 7/4 + 2c.

**3.3 Theorem C: sound.** The exact sympy steps re-run cleanly. The IBP boundary terms vanish:
- V(1) = 0 against log²g ~ log²(t−1);
- log g(t₁) = 0.

The x-variable letters are 0 < τ < s < 1 < 1/s < 1/τ, since tan(θ/2) < sin θ. My quadrature of the typed W₃ table
confirms the split to 5.7e-47 at 13 angles.

The W₂ evaluation in the solver's script uses three dilogarithm antiderivatives. Those are checked only numerically
(1e-121), so here is a **short independent proof of W₂ = ½ζ(2) − ½Li₂(c²)** (`v_proofs_sympy.py`):
1. Because log g(t₁) = 0, the moving endpoint contributes nothing.
2. ∂_c log g = 4c/(t − cos α). Hence
   dW₂/dc = −(4c/(1+cos α))[log((t₁−cos α)/(1−cos α)) − log((t₁+1)/2)] = (2/c)·log s = log(1−c²)/c,
   which equals d/dc[½ζ(2) − ½Li₂(c²)].
3. At c = 0: W₂ = −∫₀¹log v dv/(1+v) = π²/12.

sympy confirms each step. It leaves the final simplification as 2log s + 2log(1+s) − 2log(s²+s), which is 0 for s > 0.
Numerically the derivative agrees to 1e-165 at 4 points, and W₂(0) − π²/12 = 1.2e-41.

**3.4 Theorem B: labels are honest; B1 can be upgraded to PROVED.**
- The solver labels B1–B3 as EXACT-k plus NUMERICAL, and says the all-n summation was not done. That is accurate.
- **Upgrade (new, all n).** Let A_m = Σ_{k≤m}(−1)^{k+1}/k. Then A_{2n−1} = H_{2n−1} − H_{n−1} = D_n.
  - The c^N coefficient of the HPL form of B1 is −X_N + X_{N−1} + Y_N + Y_{N−1}, with X_N = 1/N² + A_{N−1}/N and
    Y_N = (−1)^{N+1}/N² − B_{N−1}/N.
  - Using B_{2n−1} = A_{2n−1}, B_{2n−2} = −A_{2n−2} and B_{2n} = −A_{2n}, this equals exactly D_n/(n(2n−1)) − 1/(2n²)
    for N = 2n, and 0 for odd N (`v_proofs_sympy.py`: both simplify to 0; the sign rules are checked for n < 60).
  - Together with the PROVED formula for μ₁,ₙ, **B1 (HPL form) holds for all α**.
  - The classical form of B1 follows exactly from H_{±1,∓1}(x) = Li₂((1∓x)/2) − Li₂(½) − log2·log(1∓x)
    (differentiate) and 2Li₂(½) = ζ(2) − log²2.
- B2 and B3 remain EXACT-600 / EXACT-640 plus NUMERICAL. A proof route for B2 would be the same nested-sum reduction:
  write the B2 coefficients and the μ₂,ₙ formula as alternating double harmonic sums and reduce by stuffle.

**3.5 The "1/c is needed" claim.**
- The task summary says 1/c-terms "are needed". result.md §3 says only "in these natural forms", and Open question 4
  asks whether some basis removes them.
- Given B2 as an identity of functions, the stronger claim is true. Iterated integrals are linearly independent over
  the rational-function field ℂ(c) (Deneufchâtel–Duchamp–Minh–Solomon 2011). So the coefficient of H_{1,0,1}(c) in M₂
  is intrinsic, and from B2 it is (1−c²)/(2c).
- Any polylogarithm of an argument in ℚ(c) expands into iterated integrals with *constant* coefficients. So a
  representation with all prefactors in ℂ ⊕ ℂc would give that coefficient in ℂ ⊕ ℂc, which is a contradiction.
- Hence the "ℚ ⊕ ℚc" rule genuinely fails at weight 3, conditional only on B2. This also answers Open question 4
  (for polylogarithmic bases): no.

**3.6 Minor wording.** "M_k ∋ (−1)^k[(1−c)Li_{k+1}(c) + (1+c)Li_{k+1}(−c)]" refers to the HPL coefficients of
H_{0…0,±1}. In the classical form B2′ the Li₃(±c) coefficients come out as (1+c) and (1−c), swapped, because other
HPLs also feed Li₃(±c) when converted. Worth one clause.

**3.7 Small items.**
- The PSLQ relations in `flat_expansion.py` (coefficients up to 5.9·10⁶) would not pass the small-coefficient standard
  on their own. They are only corroboration; the exact computation is what the claim rests on. Fine as labelled.
- The k = 4, 5 single-ζ patterns (Open question 2) are exact only for n ≤ 12, as the solver says. My re-run of
  `flat_single_zeta_k5.py` output agrees: (n+3)/(2n+3) and −3/(2n+3) after multiplying by (4n²−1). A symbolic-n re-run
  at k = 5 was started but had not finished when this report was written, because the machine was heavily loaded. No
  claim depends on it.

## 4. Literature (light check)
- One web search for large-p superellipse/Lamé perimeter expansions found only Rodríguez–Montilla, arXiv 2607.09048
  (hypergeometric series; the limit 4(a+b)). There is no 1/p corner expansion and no polylogarithmic coefficients.
- This is consistent with the solver's hedged "to my knowledge new" and with the project's earlier prior-art logs.
- The tools (HPLs, Beta derivatives) are standard, as stated.
- The open status of the ℚ-independence of 1, ζ(2), ζ(3) was confirmed (arXiv 1401.5322; Waldschmidt's survey).
- No novelty claim is contradicted, but the check was not exhaustive.

## 5. Corrections to apply to result.md
1. §2: replace "standard independence of {1, π², ζ(3)}" with "the (open) conjecture that 1, π², ζ(3) are ℚ-linearly
   independent; formally, the ζ(3)-part is canonical by the polygamma reduction".
2. §2 derivation step 2: add the complex-u (or dominated-derivative) uniformity sentence.
3. Headline: keep "as a power series in c² (flat-angle normalisation)" wherever "canonical ζ(3)-coefficient of M₂ is
   0" appears.
4. §3: B1 can be relabelled PROVED (proof in §3.4 here). W₂ has the short proof in §3.3.
5. §3/§7 Q4: the "1/c needed" statement is provable given B2 (§3.5). Change "in these natural forms" accordingly, and
   answer Open question 4 for polylogarithmic bases.
6. Clarify the leading-letter pattern as an HPL-form statement (§3.6).

## 6. Status by claim

| Claim | Verifier status |
|---|---|
| (a) M₂ at 9 angles, ~318 digits; π/2 = squircle; 3π/4, π/6 = D10b | CONFIRMED (independent method, 46 digits; reference files, 1e-322) |
| Theorem A: μ_{k,n} ∈ ℚ[ζ(2..k)]; ζ(2)/(4n²−1) in μ₂; −ζ(3)/(4n²−1) in μ₃ | CONFIRMED, PROVED (minor uniformity gap; M₃ uniqueness conditional on an open conjecture) |
| Canonical ζ(3) of M₂ = 0; of M₃ = −Φ(c) | CONFIRMED as a c²-Taylor statement (Correction 2) |
| Theorem C: complex letters cancel; W₂ = ½ζ(2) − ½Li₂(c²) | CONFIRMED, PROVED (independent W₂ proof added) |
| B1 (M₁, all α) | CONFIRMED; **upgraded to PROVED** |
| B2 / B2′ (M₂, all α), Z(c) = 7/4 + 2c, field ℚ(c) | CONFIRMED NUMERICAL + EXACT to n ≤ 300 (my code); not proved for all n |
| B3 (M₃ closed form, Ψ) | CONFIRMED NUMERICAL (5 fresh angles, 46 digits); Ψ exact for n ≤ 40 (my code) |
| 2π/3 closed form; bases with ζ(3)-coefficient 3/16, 123/16, −29/4 | CONFIRMED (320 digits; identities checked) |
| π/3 table, ζ(3)-coefficient 7/4 + √3 (uniform basis) | CONFIRMED (B2′ at c = √3/2 to 4.8e-322; arguments and coefficients re-derived by hand) |
| π/2 minimal-basis coefficient 189√2/64 − 19/16 | CONFIRMED (re-run; the note's identities checked at 60 digits in the script) |
| Novelty | not contradicted (light check) |

## 7. Verifier scripts (in `verify/`)
- `v_cauchy.py`: independent M_k by Cauchy contour; outputs `vc_*.txt` (run A) and `runB_vc_*.txt` (run B).
- `v_compare_values.py`: run A vs run B vs solver values vs squircle.
- `v_flat_numeric.py`: μ_{k,n} from the original g-integral, compared with the exact table, plus PSLQ.
- `v_flat_symbolic.py`: own symbolic ζ-coefficients and the Ψ check.
- `v_taylor_exact.py`: exact Taylor check of B1 and B2, n ≤ 300.
- `v_closed_forms.py`: hand-typed B0, B1, B2′, HPL B2, W₂, W₃, 2π/3 and π/3, at 13 angles and at 330 digits.
- `v_m3_newangles.py`: B3 at fresh angles.
- `v_flat_sum.py`: partial sums vs independent values.
- `v_basis_change.py`: the 2π/3 bases.
- `v_d10b_compare.py`: comparison with D10b.
- `v_proofs_sympy.py`: the W₂ and B1 proofs.
- Each script has a `.out` file. Scripts read the solver's data files (`values_*.txt`, `flat_exact_N*_K3.json`,
  `hpl_fit_*.json`) from the scratch copy, so run them from the scratch folder or alongside the solver's files.

# G1: ζ(3) at a corner of any angle (the p → ∞ end)

Folder: `research/10_frontier/G1_zeta3_any_angle/`. Python: `.venv/bin/python`. Every number below is printed by a
script in this folder (the `.out` file is named after the script).

**Notation.** Interior angle α ∈ (0, π), c = cos(α/2), s = sin(α/2), τ = tan(α/4).
The corner function (as in `corner_series.py`) is
M(u;α) = ∫₀¹ [(1+g)^(−1−u) √(1+g^(2−2u) − 2cos α·g^(1−u)) − 1] g^(u−1) dg = Σ_k M_k(α) u^k.
H_w(c) are harmonic polylogarithms (HPLs; Remiddi–Vermaseren) with letters 0, 1, −1, so that H₁ = −log(1−c),
H₋₁ = log(1+c), H_{0,0,1} = Li₃(c), H_{a,w}(c) = ∫₀^c f_a H_w, with f₀ = 1/t, f₁ = 1/(1−t), f₋₁ = 1/(1+t).
I write Φ(c) = Σ_{n≥1} c^{2n}/(4n²−1) = ½ − ((1−c²)/(2c))·artanh c.

Labels used below:
- **PROVED**: exact sympy algebra, or a complete argument written out here.
- **EXACT-k**: the first k Taylor coefficients agree in exact rational arithmetic.
- **NUMERICAL**: high-precision agreement, with digit counts given.
- **CONJECTURE**: everything else.

---

## 0. Short answers

**(a)** M₂(α) is computed to about 318 digits at 9 angles: π/2, 2π/3, π/3, 3π/4, π/6, 5π/6, π/4, and the generic
angles 1.3 and 2.0 rad (table in §5).
- π/2 matches the squircle note's 1200-digit value to 3.3·10⁻³²² (`squircle_check.out`).
- 3π/4 and π/6 match D10b's 330-digit values to 10⁻³²².

**(b) Closed forms.**
- I found a closed form for **every** angle (Theorem B2). It contains Li₃ at 8 arguments, each a rational function of c.
- 2π/3 (over ℚ, §3) needs only Li₃(1/3), Li₃(−1/2) and ζ(3).
- π/3 (over ℚ(√3), §3) is the specialisation c = √3/2 of Theorem B2.

**(c) The coefficient of ζ(3).**
1. **Canonically it is 0 for every α (PROVED, Theorem A).** Expand about the flat angle α = π, where M ≡ 0. Then
   M₂(α) = Σ_{n≥1} μ_{2,n} c^{2n}, with every μ_{2,n} ∈ ℚ + ℚζ(2), and with the explicit formula
   μ_{2,n} = r_n + ζ(2)/(4n²−1). No ζ(3) and no log 2 occur in any Taylor coefficient. So the intrinsic
   decomposition is M₂ = ζ(2)Φ(c) + R₂(c), where R₂ has rational Taylor coefficients.
   - The ζ(3) seen in closed forms comes from evaluating these trilogarithmic functions of c at a particular c.
   - Its coefficient in such a closed form depends on the basis. At 2π/3 the same number has ζ(3)-coefficient
     3/16, 123/16 or −29/4 in three equally valid bases (`m2_2pi3_pslq.out`).
2. **In the uniform classical basis it is the polynomial Z(c) = 7/4 + 2c (NUMERICAL, 300+ digits at 9 angles).**
   This is the basis of Theorem B2: Li₃ at the 8 standard HPL arguments
   c, −c, 1−c, 1/(1+c), (1+c)/2, (1−c)/(1+c), −(1−c)/(1+c), 2c/(1+c), all in (−1,1).
   - Z is forced by the Li₃ coefficients: Z(c) = −Σ_i C_i(c)·[ζ(3)-part of Li₃(a_i(0))]. This is the flat
     normalisation, checked symbolically (`zeta3_bookkeeping.out`).
   - Reducing B2 at c = 1/√2 with the note's Lemma (lem:trilog) reproduces the squircle coefficient
     189√2/64 − 19/16 exactly.
   - Reducing B2 at c = 1/2 to the basis {Li₃(1/3), Li₃(−1/2)} gives 3/16, which agrees with an independent PSLQ.
3. **Where ζ(3) is intrinsic: in M₃, with coefficient −Φ(c) (PROVED, Theorem A).** The coefficient of c^{2n} in M₃
   is r′_n + q′_nζ(2) − ζ(3)/(4n²−1). The same elementary function Φ multiplies +ζ(2) in M₂ and −ζ(3) in M₃.
   - At the square corner the canonical ζ(3) content of M₃ is −ζ(3)(½ − (√2/4)log(1+√2)).

---

## 1. Theorem C: exact weight split; the complex letters cancel (PROVED)

Script: `m2_structure.py` (output `m2_structure.out`), `m2_w2_identity.py`, `m2_split_check.py`.

1. With the Euler substitution t = g + √(1−2cos α·g+g²) (D10c), the integrand becomes
   M₂ = ∫₁^{t₁} [R_A log²g + R_B log g·log(1+g) + R_C log²(1+g)] dt, with t₁ = 1+2s.
   - This was re-derived independently: the u-expansion treats L = log g and m = log(1+g) as symbols.
   - R_B and R_C agree with D10c.
2. **Key identity (sympy, exact):** R_A = R_A^s − s²·d/dt[(t+1)/P], with P = t² − 2cos α·t + 1. Here R_A^s has
   only simple poles, at the real letters 1, −1, cos α, r₁,₂ = −1 ± 2c. Its residues are
   (0, −1, ½, −(c²+1)/(4c), +(c²+1)/(4c)).
3. **Integration by parts.** Take V = ½ − s²(t+1)/P = (t−1)(t−cos α)/(2P). Then V(1) = 0 and log g(t₁) = 0, and
   V·(d/dt)log g = 1/(2(t+1)) **exactly**. Hence ∫ R_A^c log²g dt = −∫ log g dt/(t+1).
4. **Consequence.** The double poles at e^{±iα} that D10c found (its open question 4) contribute only a weight-2 term:
   **M₂(α) = W₂ + W₃**, where
   - **W₂ = −∫₁^{t₁} log g dt/(t+1) = ½ζ(2) − ½Li₂(c²)** (PROVED: the log bookkeeping is done by sympy in
     `m2_w2_identity.py`, then Landen and Euler reflection). At π/2 this is the note's (π²+6log²2)/24.
   - **W₃** is a pure weight-3 iterated integral on the real alphabet only. In x = (t−1)/(2s) ∈ [0,1]:
     W₃ = ∫₀¹ Σ_σ [A_σ ℓ_g² + B_σ ℓ_g ℓ_m + C_σ ℓ_m²] dx/(x+λ_σ),
     ℓ_g = log(x/s) + log(1+sx) − log(1+x/s),  ℓ_m = log(1+x/τ) + log(1+τx) − log(1+x/s).
   - The poles and residues are:

     | λ_σ | (A_σ, B_σ, C_σ) |
     |---|---|
     | 0 | (0, −1, ½) |
     | 1/s | (−1, 1, −½) |
     | s | (½, 0, ½) |
     | τ | (−(c²+1)/(4c), c, −c) |
     | 1/τ | (+(c²+1)/(4c), −c, c) |

     The ordering 0 < τ < s < 1 < 1/s < 1/τ holds for every α. So there are no branch case distinctions.
5. Check: W₂ + W₃ by quadrature equals the series M₂ to 150 digits at all 9 angles (`m2_split_check.out`).

---

## 2. Theorem A: expansion about the flat angle; the canonical zeta content (PROVED)

Scripts: `flat_exact.py` (exact μ_{k,n}, n ≤ 320), `flat_symbolic_proof.py` (symbolic n), `flat_mu2_symbolic.py`,
`flat_expansion.py` (independent PSLQ), `flat_sum_check.py`.

**Derivation (complete).**
1. Write 1 + h² − 2cos α·h = (1+h)² − 4c²h, with h = g^(1−u). Then
   (1+h)√(1−4c²h/(1+h)²) = (1+h) − Σ_{n≥1} κ_n 4ⁿ c^{2n} hⁿ(1+h)^(1−2n), where κ_n4ⁿ = C(2n,n)/(2n−1).
2. **M(u;π) = 0 identically**, because ∫₀¹(1+g)^(−1−u)g^(u−1)dg = 2^(−u)/u. Hence
   M(u;α) = −Σ_n κ_n4ⁿ c^{2n} J_n(u). The series converges uniformly for c² < 1, since 4c²h/(1+h)² ≤ c².
3. Each term is invariant under g ↦ 1/g, so J_n = ½∫₀^∞. Substituting y = g/(1+g) gives
   J_n(u) = ½∫₀¹ y^((n−1)(1−u)) (1−y)^(n−1+nu) [1 + y(((1−y)/y)^u − 1)]^(1−2n) dy.
4. The u^k coefficient is a polynomial in y, log y and log(1−y). So every μ_{k,n} is a finite sum of derivatives
   ∂_p^i∂_q^l B(p,q) at the integers p = n+j, q = n. These are Bell polynomials in the cumulants of log B, that is
   in polygammas at integers: ψ(N)+γ = H_{N−1}, and ψ^(m)(N) = (−1)^(m+1)m!(ζ(m+1) − H^(m+1)_{N−1}).
5. **Hence μ_{k,n} ∈ ℚ[ζ(2), …, ζ(k)], of weight ≤ k, for all k and n.**
   - The pure cumulants κ_{r,0} and κ_{0,s} are rational, because the ζ's cancel.
   - Only the mixed cumulants κ_{r,s} = −ψ^(r+s−1)(p+q), with r, s ≥ 1, carry ζ(r+s).

**Results (symbolic n, `flat_symbolic_proof.out`):**
- **μ_{2,n} = r_n + ζ(2)/(4n²−1), with r_n ∈ ℚ.** No ζ(3) can occur, because [u²] has log-degree 2. Explicitly
  (`flat_mu2_symbolic.out`, which checks it against the exact table for all n ≤ 320):

  μ_{2,n} = ζ(2)/(4n²−1) − D_n²/(2n(2n−1)) + D_n/(2n²) + H^(2)_{n−1}/(2n(4n²−1)) − H^(2)_{2n−1}/(2n(2n−1))
  − (4n²−2n−1)/(4n³(2n+1)),

  where D_n = H_{2n−1} − H_{n−1}.
- μ_{0,n} = −1/(n(2n−1)), and μ_{1,n} = D_n/(n(2n−1)) − 1/(2n²).
- **The coefficient of ζ(3) in μ_{3,n} is −1/(4n²−1)**, as a symbolic identity in n.

**Corollary (canonical decomposition; PROVED).**
- **M₂(α) = ζ(2)·Φ(c) + R₂(c)**, with R₂ ∈ ℚ[[c²]].
- **M₃(α) = −ζ(3)·Φ(c) + ζ(2)·Ψ(c) + R₃(c)**, with Ψ, R₃ ∈ ℚ[[c²]].

The splitting of each Taylor coefficient into rational and ζ-parts is unique, given the irrationality of π² (and,
for M₃, the standard independence of {1, π², ζ(3)}). In this canonical sense the ζ(3)-coefficient of M₂ is
identically **0**.

**Checks.**
- An independent PSLQ of μ_{k,n} (n ≤ 3, 120 digits, against the full alternating-MZV basis of weight ≤ k+1)
  confirms that no log 2, log²2, π²log 2 or log³2 occurs (`flat_expansion_test.out`).
- Σ_{n≤320} μ_{k,n}c^{2n} reproduces the independent series values to the truncation error (`flat_sum_check.out`):
  10⁻³²² at 5π/6, 10⁻²⁷⁴ at 3π/4, 10⁻¹⁹⁹ at 2π/3, 10⁻¹⁰² at π/2.

---

## 3. Theorem B: general-angle closed forms

**Method.** The exact Taylor coefficients (§2) are fitted, in exact rational arithmetic, against
{c^e H_w(c)}, with w a word over {0,±1} with nonzero last letter and e ∈ {−1,0,1} (for M₃, e ∈ {−2,…,1}).
- The systems are heavily overdetermined: M₂ has 80 unknowns and 401 equations; M₃ has 322 unknowns and
  641 equations.
- A solution was accepted only if it satisfied every equation (`hpl_fit.py`).
- The resulting identities are **EXACT through c^400 (M₂, n ≤ 200) and c^640 (M₃)**, and **NUMERICAL** beyond that.
- An all-n proof would need the standard summation of the six harmonic-sum generating functions in the μ_{2,n}
  formula above, or a holonomic certificate. **I did not carry this out.**

**B0 (classical, known).** M₀ = (1−c)H₁ − (1+c)H₋₁ = −(1+c)log(1+c) − (1−c)log(1−c).

**B1 (new compact M₁, all α).**
M₁(α) = ζ(2) − log²2 − log 2·M₀(α) − Σ_{ε=±1}(1−εc)[Li₂(εc) + Li₂((1−εc)/2)].
- HPL form: M₁ = −(1−c)[H_{0,1} + H_{1,−1}] + (1+c)[H_{0,−1} − H_{−1,1}].
- It replaces D10c's 31-Li₂ formula, and its coefficients lie in ℚ ⊕ ℚc.
- Checks: 320 digits at 9 angles (`m1_classical_check.out`, `report_formula_check.out`).

**B2 (M₂, all α; HPL form).**

M₂(α) = ½ζ(2) − ½Li₂(c²) − ζ(2)·((1−c²)/(2c))·artanh c + S(c) + S(−c),

S(c) = (1−c)[H_{0,0,1} + H_{0,1,−1} − H_{1,−1,1}](c) + ((1−c²)/(2c))H_{1,0,1}(c) − ((1−c)²/(2c))H_{1,0,−1}(c).

Here S(−c) means the same expression with every H_w(−c) replaced by (−1)^{#nonzero letters} H_{−w}(c). It is
written out in `report_formula_check.py`. **There is no ζ(3) term.**
- Checks: 300 digits at all 9 angles (`hpl_eval_300.out`, residuals ≤ 3·10⁻³¹⁶). This includes c = 0.966 (π/6),
  where the Taylor tail beyond the fitted range is about 10⁻⁶. So the check probes coefficients the fit never saw.

**B2′ (M₂, all α; uniform classical form, real-valued with principal branches).** Write ℓ₊ = log(1+c),
ℓ₋ = log(1−c), ℓ₀ = log c, λ = log 2 and q = (1−c²)/(2c). Then

M₂ = ζ(2)Φ(c) − Li₂(c) − Li₂(−c)
 + (1+c)Li₃(c) + (1−c)Li₃(−c) + ((1−c)/c)Li₃(1−c) − ((1+c)/c)Li₃(1/(1+c)) − 2c·Li₃((1+c)/2)
 + (1−c)[Li₃((1−c)/(1+c)) + Li₃(−(1−c)/(1+c))] − 2c·Li₃(2c/(1+c))
 + Li₂((1+c)/2)[(1+c)ℓ₊ − (1−c)ℓ₋] + q[Li₂(c) − Li₂(−c)](ℓ₊+ℓ₋)
 − (1−c)ℓ₋²ℓ₊ + (1−c)ℓ₋ℓ₊² + ((4c²−c+1)/(6c))ℓ₊³ + qℓ₀ℓ₋² − (1−c)ℓ₀ℓ₋ℓ₊ − ((1+c²)/(2c))ℓ₀ℓ₊²
 + λ[(1−c)ℓ₋² − (1−c)ℓ₋ℓ₊ − cℓ₊²] + λ²((1−c)/2)(ℓ₊−ℓ₋) − π²[((1−c)/(6c))ℓ₋ + ((1+c)/(6c))ℓ₊]
 **+ (7/4 + 2c)ζ(3)** − (c/6)π²λ + (c/3)λ³.

Checks: 320 digits at 9 angles (`m2_classical_assemble.out`), and the typed formula is re-checked at 120 digits
(`report_formula_check.out`).
- The 10 HPL→classical identities used here (`classical_dict.out`) are standard. Each was re-derived by a
  function-level fit with a 40-function basis of full rank, then checked at 40 fresh points in (0, 0.55) to
  10⁻¹¹⁰, and at x = 0.7, 0.9, 0.97 to 10⁻⁴⁵.
- Coefficient ring: the coefficients lie in c⁻¹(ℚ + ℚc + ℚc²). **The D10b/D10c rule "ℚ ⊕ ℚc" does not extend to
  weight 3** in these natural forms: the 1/c-terms (1±c)/c, (1−c²)/(2c), (4c²−c+1)/(6c) are needed. This matches
  the R_A residues ∓(c²+1)/(4c) of §1. The field rule K = ℚ(c) does hold: all arguments lie in ℚ(c), and s never
  appears.

**B3 (M₃, all α; bonus).**
M₃ = −ζ(3)Φ(c) + ζ(2)Ψ(c) + R₃(c), where
- Ψ = ½ + (¾ − 1/(4c) − 1/(2c²))H₁ + (−¾ − 1/(4c) + 1/(2c²))H₋₁ + ((1−c²)/(4c))[H_{0,1} + H_{0,−1} + H_{1,−1} − H_{−1,1}];
- R₃ has 74 HPL terms of weights 2–4 (`hpl_fit_k3_rational.out`).

Status: EXACT through c^640; NUMERICAL at 6 angles to 58–61 digits at 60-digit working precision
(`m3_hpl_check_*.out`). This includes c = 0.966, where the tail beyond the fit is about 2·10⁻¹⁰, and it reproduces
the squircle M₃ = 0.02277052149178657745673…

The leading-letter pattern is M_k ∋ (−1)^k[(1−c)Li_{k+1}(c) + (1+c)Li_{k+1}(−c)] for k = 0, 1, 2, 3.

**Special angles.**
- **2π/3 (c = ½), over ℚ, minimal basis** (`m2_2pi3_pslq.out`). PSLQ was run on W₃ only, at 300 digits, against
  the 11-element basis of real weight-3 periods of MT(ℤ[1/6]). The relation has maximum coefficient 408 and
  residual 2·10⁻³²¹. The formula is re-checked at 120 digits in `report_formula_check.out`.

  M₂(2π/3) = Li₂(⅓) + log²2 − log2·log3 + ½log²3
   + (15/2)Li₃(⅓) + (17/2)Li₃(−½) + **(3/16)ζ(3)** − (5/12)π²log2 + ¼π²log3 − (1/12)log³2 + (5/4)log²2·log3
   − (7/4)log2·log²3 + (5/8)log³3 − (7/2)Li₂(⅓)log2 + (15/4)Li₂(⅓)log3
   = −0.009567895652956384775719722385651527995…

  The same number in other bases:

  | Basis | ζ(3)-coefficient |
  |---|---|
  | {Li₃(2/3), Li₃(−½)} | 123/16 |
  | {Li₃(⅓), Li₃(¼)} | −29/4 |

- **π/3 (c = √3/2), over ℚ(√3).** This is B2′ at c = √3/2.

  | Li₃ argument | coefficient |
  |---|---|
  | √3/2 | 1 + √3/2 |
  | −√3/2 | 1 − √3/2 |
  | 1 − √3/2 | 2√3/3 − 1 |
  | 4 − 2√3 | −1 − 2√3/3 |
  | (2+√3)/4 | −√3 |
  | 7 − 4√3 | 1 − √3/2 |
  | −(7−4√3) | 1 − √3/2 |
  | 4√3 − 6 | −√3 |

  - ζ(3) coefficient: 7/4 + √3.
  - Weight-2 part: ½ζ(2) − ½Li₂(¾).
  - The logs are of (2±√3)/2 and √3/2 (full list in `m2_classical_assemble.out`). Value −0.04287175497685288253…
  - A reduction to a minimal ℚ(√3) basis (S = {2,3}) was not attempted.
- **π/2.** B2′ reduces to the note's Theorem M₂, including its ζ(3)-coefficient 189√2/64 − 19/16
  (`zeta3_bookkeeping.out`).

---

## 4. The ζ(3)-coefficient, angle by angle

"Canonical" means the flat-corner and HPL normalisation of Theorem A. "Uniform" means basis B2′, where
Z = 7/4 + 2c. "Minimal" means a small field-specific basis.

| α | c | canonical | uniform Z(c) | minimal basis |
|---|---|---|---|---|
| π/2 | 1/√2 | 0 | 7/4 + √2 = 3.16421… | 189√2/64 − 19/16 (note's basis {Li₃(b), Li₃(−b), Li₃(1/√2)}) |
| 2π/3 | ½ | 0 | 11/4 | 3/16 (basis A {Li₃(⅓), Li₃(−½)}); 123/16 (B); −29/4 (C) |
| π/3 | √3/2 | 0 | 7/4 + √3 | — |
| 3π/4 | cos(3π/8) | 0 | 7/4 + √(2−√2) | — |
| π/6 | cos(π/12) | 0 | 7/4 + (√6+√2)/2 | — |
| 5π/6, π/4, 1.3, 2.0 | | 0 | 7/4 + 2c (values in `summary_table.out`) | — |

The canonical coefficient of ζ(3) in M₃ (PROVED) is −Φ(c):

| α | −Φ(c) |
|---|---|
| π/2 | −½ + (√2/4)log(1+√2) = −0.18839… |
| 2π/3 | −½ + (3/8)log 3 = −0.08802… |
| π/3 | −½ + log(2+√3)/(4√3) = −0.30991… |

**Is it a simple function of c?**
- Canonically it is 0 for M₂, which is trivially polynomial.
- In the uniform classical basis it is the linear polynomial 7/4 + 2c.
- Other uniform choices of Li₃ representatives shift it by ℤ-combinations of the Li₃ coefficients C_i(c), which
  are Laurent polynomials in ℚ[c, 1/c]. For example, trading Li₃(1/(1+c)) for Li₃(c/(1+c)) adds −(1+c)/c.
- So no choice-free polynomial "ζ(3)-coefficient" exists beyond the canonical one.

---

## 5. Evidence table

| Check | Label | Result |
|---|---|---|
| `m2_numerics.py`: M₀…M₃ at 9 angles, dps 320 vs 350 | NUMERICAL | diffs ≤ 2·10⁻³²⁰ (M₂), ≤ 4·10⁻³¹⁸ (M₃) |
| M₂(π/2) vs squircle 1200-digit value | NUMERICAL | 3.3·10⁻³²² (`squircle_check.out`) |
| M₂(3π/4), M₂(π/6) vs D10b 330-digit values | NUMERICAL | 1.4·10⁻³²³, 2.0·10⁻³²² |
| `m2_structure.py`: integrand, key identity, IBP identity, residues | PROVED | all simplifications = 0 |
| `m2_w2_identity.py`: W₂ = ½ζ(2) − ½Li₂(c²) | PROVED | log part = log s·log c²; numerical 10⁻¹²¹ at 9 angles |
| `m2_split_check.py`: W₂ + W₃ (quadrature) vs series | NUMERICAL | ≤ 3·10⁻¹⁵¹, 9 angles |
| `flat_symbolic_proof.py`: ζ(2)-coefficient of μ_{2,n} and ζ(3)-coefficient of μ_{3,n} | PROVED | 1/(4n²−1) and −1/(4n²−1) |
| `flat_exact.py`: exact μ_{k,n}, k ≤ 3, n ≤ 320 | PROVED (exact) | same pattern for every n ≤ 320 |
| `flat_expansion.py`: PSLQ of μ_{k,n} in the full alt-MZV basis | NUMERICAL | max coefficient ≤ 5.9·10⁶ for n ≤ 3, residual ≤ 2·10⁻¹³⁶; no log 2 or weight-3 constants in μ₂ |
| `flat_sum_check.py`: Σ_{n≤320} μ c^{2n} vs series | NUMERICAL | down to truncation error (10⁻³²² … 10⁻¹⁰²) |
| `flat_mu2_symbolic.py`: harmonic-number formula for μ_{1,n}, μ_{2,n} | PROVED | matches the exact table for n ≤ 320 |
| `hpl_fit.py`: M₀, M₁, Φ, R₂ (n ≤ 200); M₃ parts (n ≤ 320) | EXACT-400 / EXACT-640 | rank(A) = rank([A\|b]) = #unknowns |
| `hpl_eval_par.py`: B0, B1, B2 (HPL quadrature) vs series | NUMERICAL | ≤ 4·10⁻³¹⁶, 9 angles |
| `m2_classical_assemble.py`: B2′ vs series | NUMERICAL | ≤ 1.3·10⁻³²⁰, 9 angles |
| `m1_classical_check.py`: B1 | NUMERICAL | ≤ 5·10⁻³²¹ |
| `classical_dict.py`: 10 HPL→Li₃ identities | NUMERICAL | 10⁻¹¹⁰ (40 fresh points), 10⁻⁴⁵ near 1 |
| `m2_2pi3_pslq.py`: M₂(2π/3), 12-element basis, 300 digits | NUMERICAL | max coefficient 408, residual 2·10⁻³²¹ (3 bases) |
| `zeta3_bookkeeping.py`: flat normalisation of Z; squircle and 2π/3 reductions | PROVED (symbolic) / NUMERICAL | 0, 0, 0; 189√2/64 − 19/16; 3/16 |
| `report_formula_check.py`: formulas as typed in this file | NUMERICAL | 10⁻¹²¹ (B1, B2′, 2π/3), 10⁻⁶¹ (HPL B2) |
| `m3_hpl_check.py`: B3 vs series at 6 angles | NUMERICAL | 5·10⁻⁶¹ … 1.3·10⁻⁵⁸ (60-digit working precision) |
| `flat_single_zeta_k5.py`: single-ζ(k) coefficient of μ_{k,n}, k ≤ 5, n ≤ 12 | PROVED (exact) for these n | k=2: 1/(4n²−1); k=3: −1/(4n²−1); k=4: (n+3)/((2n+3)(4n²−1)) (formal; ζ(2)² is also present); k=5: −3/((2n+3)(4n²−1)) |

**Values** (318 digits in `values_*.txt`; table in `summary_table.out`):

| α | M₂(α) |
|---|---|
| π/2 | −0.02216191108726939942375493773979151497795504632874345620787 |
| 2π/3 | −0.00956789565295638477571972238565152799500830163335114401714 |
| π/3 | −0.04287175497685288253230642505250886346370973319237454663412 |
| 3π/4 | −0.00537310154272117376549989171498160264618991637117763540589 |
| π/6 | −0.07614050794643136197265967844532319457760068500736669055903 |
| 5π/6 | −0.00239241339321074168609894060876908809283689932634121310568 |
| π/4 | −0.05775219575871169616220844109169304681945137175477034593295 |
| 1.3 rad | −0.03157868954628429226452796746009672029619277608168643213317 |
| 2.0 rad | −0.01139513052622226513790857911222593406509271519964252237458 |

---

## 6. Literature

- A web search, plus the project's earlier prior-art logs (`05_final_verify/corner_angle/*queries.log`,
  `kimi_answer.md`), found no treatment of the p-smoothed corner coefficients M_k(α). The searches covered
  superellipse and Lamé-curve perimeter expansions, and amoeba or limit-shape boundary lengths.
- Rodríguez–Montilla (arXiv 2607.09048) give hypergeometric series for the Lamé perimeter and the p → ∞ limit 4(a+b),
  without the 1/p corner coefficients.
- The tools are standard: HPLs (Remiddi–Vermaseren, IJMPA 15 (2000) 725), and Beta-function derivatives and
  polygammas at integers.
- To my knowledge the following are new, but this is not exhaustively checked:
  - the flat-angle expansion (Theorem A);
  - the canonical pairing +ζ(2)Φ / −ζ(3)Φ;
  - the cancellation of the complex letters;
  - the closed forms B1–B3.

## 7. Open questions

1. **All-n proof of B1–B3.** Sum the explicit harmonic-number formula for μ_{2,n} (§2) into HPLs, or give a
   holonomic certificate. The Taylor identity is exact to c^400 and numerically confirmed far beyond.
2. **The pattern in general. The naive guess fails.** The single-ζ(k) coefficient of μ_{k,n} is ±1/(4n²−1) for
   k = 2, 3. For k = 4 and 5 it is (n+3)/((2n+3)(4n²−1)) and −3/((2n+3)(4n²−1)) (`flat_single_zeta_k5.out`; exact for
   n ≤ 12, pattern guessed). So M₅'s canonical ζ(5)-coefficient is presumably −3Σc^{2n}/((2n+3)(4n²−1)), again
   elementary in c. Open: prove the general-k formula, and explain why the same Φ carries +ζ(2) in M₂ and −ζ(3) in M₃.
3. **Minimal-basis reduction of M₂(π/3) over ℚ(√3)** (S = {2,3}; 8 Lyndon trilogs plus 2 primitive weight-3 periods).
   Do the same for the quartic angles (3π/4, π/6), where the ζ(3)-coefficient in a minimal basis is not yet known.
4. **Interpretation of Z(c) = 7/4 + 2c.** It equals −Σ C_i(c)·z(a_i(0)), so it is fixed by the Li₃ coefficients
   and the arguments' values at the flat angle. Is there a natural basis making the 1/c coefficients disappear?
5. **The other end.** The flat-angle expansion has radius |c| < 1; the cusp α → 0 is c → 1. Find the asymptotics
   of M_k as α → 0, and whether the constant term there carries a canonical ζ(3).

## 8. Scripts

| Script | Purpose |
|---|---|
| `corner_series.py` | copy of the reference series code |
| `m2_numerics.py` | values |
| `m2_structure.py` | Theorem C, exact algebra |
| `m2_split_check.py` | W₂ + W₃ numerics |
| `m2_w2_identity.py` | W₂ = ½ζ(2) − ½Li₂(c²) |
| `flat_expansion.py` | PSLQ of the flat coefficients |
| `flat_exact.py` | exact flat coefficients, ℚ[ζ₂, ζ₃, ζ₄] |
| `flat_symbolic_proof.py` | symbolic-n proof of the ζ-coefficients |
| `flat_single_zeta_k5.py` | single-ζ coefficients for k ≤ 5 |
| `squircle_check.py` | M₂(π/2) vs the squircle value |
| `flat_mu2_symbolic.py` | harmonic-number formula |
| `flat_sum_check.py` | partial-sum check |
| `hpl_fit.py` | exact HPL fits |
| `hpl_eval.py`, `hpl_eval_par.py` | HPL numerics, 300-digit checks |
| `classical_dict.py` (+ `classical_nullspace.py`) | HPL → Li₃ dictionary |
| `m2_classical_assemble.py` | uniform classical M₂ and Z(c) |
| `m1_classical_check.py` | B1 |
| `m2_2pi3_pslq.py` | 2π/3 minimal forms |
| `zeta3_bookkeeping.py` | basis dependence of Z |
| `m3_hpl_check.py` | B3 |
| `summary_table.py` | table |
| `report_formula_check.py` | formulas as typed here |

The JSON files hold the exact fits (`hpl_fit_k*_*.json`), the classical dictionary (`classical_dict.json`) and the
exact flat coefficients (`flat_exact_N*_K3.json`).

Note: `planner_check_W2.py` was added to this folder by the planner, not by me. It independently checks W₂ and Φ at 40 digits; running it gives differences at the 1e-40 level.

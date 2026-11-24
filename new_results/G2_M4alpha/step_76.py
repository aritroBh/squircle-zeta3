cwd = "<project>/research/08_theory/G2_M4alpha"
md = '''# G2: M4 at every corner angle, in harmonic polylogarithms, and where it needs depth two

Folder: `research/08_theory/G2_M4alpha/`. Python: the project `.venv`. Every number below is printed by a
script in this folder (the `step_NN.py` files and their outputs).

**Notation** (as G1). Interior angle α ∈ (0, π), c = cos(α/2). The corner function
M(u;α) = ∫₀¹ [(1+g)^(−1−u) √(1+g^(2−2u) − 2cos α·g^(1−u)) − 1] g^(u−1) dg = Σ_k M_k(α) u^k.
H_w(c): HPLs with letters {0, ±1}; H_1 = −log(1−c), H₋₁ = log(1+c), H_{0,0,1} = Li₃(c).
Flat-corner expansion about α = π (G1 Theorem A): M_k(α) = Σ_{n≥1} μ_{k,n} c^{2n}.

Labels: **PROVED** (exact algebra; all-n identities proved by the Z-sum method of
`W1_corner_family/B_taylor_pi4/b3_allN_proof.py`), **EXACT through c^N** (exact rational agreement of the
first N Taylor coefficients), **NUMERICAL** (high-precision agreement, digit counts given).

---

## 1. Task 1: the exact coefficients μ_{4,n}

**Which constants occur (PROVED, `step_05.py`).** μ_{4,n} ∈ ℚ + ℚζ(2) + ℚζ(3) + ℚζ(4) + ℚζ(2)², and
nothing else: the flat-corner construction expresses μ_{k,n} as a ℚ(n)-combination of Bell polynomials of
degree ≤ k in the cumulants κ_m ∼ ψ^(m−1) (m ≤ k); for k = 4 the zeta monomials are exactly
{1, ζ(2), ζ(3), ζ(4), ζ(2)²}. **No ζ(5) and no ζ(2)ζ(3) occur** (they would need degree ≥ 5). Since
ζ(2)² = (5/2)ζ(4), the π⁴-part is rational; we keep the {ζ(4), ζ(2)²} split because it is the one the
cumulants produce.

- `flat_exact.py` (G1's code, unmodified) run at k = 4, N = 640: `flat_exact_N640_K4.{json,out}`.
- `step_05.py`: independent symbolic formula for μ_{4,n} (51 monomials in H^(r)_{n−1}, H^(r)_{2n−1},
  r ≤ 4, over ℚ(n); `mu4_formula.pkl`). It agrees with `flat_exact.py` **exactly for all n ≤ 640**.

**Patterns (all PROVED from the symbolic formula, `step_07.py`, `step_11.py`).**
- [ζ(4)] μ_{4,n} = (n+3)/((2n+3)(4n²−1));  [ζ(2)²] μ_{4,n} = −n/((2n+3)(4n²−1)) — exact rational-function
  identities extracted from `mu4_formula.pkl`.
- [ζ(3)] μ_{4,n} = −[ζ(2)] μ_{3,n} — hence the ζ(3)-part of B4 is −ζ(3)·Ψ(c) with **the same Ψ** as B3's
  ζ(2)-part (G1 `hpl_fit_k3_(1,0,0).json`); PROVED for all n as part of the step_11 canonical identity.
  (This extends G1's pairing ±1/(4n²−1) at the top zeta: here a full lower-k coefficient reappears.)

**Numerical check (`step_06.py`, NUMERICAL).** Σ_{n≤400} μ_{4,n} c^{2n} vs
- the squircle value M4 = 0.0060279931126091225045161082508776… (M4a 800 digits): diff **2.2·10⁻¹⁰²**
  (truncation-limited);
- `corner_series` quadrature at α = 2π/3: diff 1.5·10⁻⁹⁸; at α = 1.3: diff 1.5·10⁻⁸⁵ (both quadrature-limited).

## 2. Task 2: the HPL form B4 — PROVED for all n

**Theorem B4.** For |c| < 1,
  M₄(α) = R₄(c) + ζ(2)·Ψ₂(c) − ζ(3)·Ψ(c) + ζ(4)·F₁(c) + ζ(2)²·F₂(c),
where
- R₄ = `hpl_fit_k4_(0,0,0).json`: **109 terms**, HPL weights 2–5, prefactors c^e, e ∈ {−3,…,1}
  (weights 2: 3 terms; 3: 8; 4: 28; 5: 70, in the even-sector combination c^e(H_w ± H_{−w}));
- Ψ₂ = `hpl_fit_k4_(1,0,0).json`: **22 terms**, weights 1–3 (needs c^{−3} at weight 3);
- Ψ = B3's ζ(2)-part (15 terms), so the ζ(3)-part is **−ζ(3)·Ψ(c)**;
- F₁, F₂ elementary (`b4_pi4_part.json`), with artanh c = (H₁ + H₋₁)/2:
    F₁ = (7c/16 − 5/(8c) + 3/(16c³))·artanh c + 9/16 − 3/(16c²),
    F₂ = (−c/16 − 1/(8c) + 3/(16c³))·artanh c + 1/16 − 3/(16c²).

**Proof status.** **PROVED for all n** (`step_11.py`): for each of the five parts, the Taylor identity
[c^N](B4 − Σμ_{4,n}c^{2n}) = 0 holds as an exact identity of ℚ(n)-combinations of alternating Z-sums,
separately for even N (equals the part of μ_{4,n}) and odd N (equals 0), with the monomial coefficients
c^{−3}…c³ checked directly. Both sides converge on |c| < 1, so B4 is an identity of analytic functions there
(and continuous up to c → 1⁻). The exponent range e ∈ {−3,…,1} is exactly the one forced by the
denominators (2n±1, 2n+3, n, n+1) of the μ-formula.
- Independently: **EXACT through c¹²⁸⁰** (`step_10.py`): all five parts' combinations reproduce the exact
  target coefficients through c^1280 (641 even equations), all odd coefficients vanish exactly; the
  rational and ζ(2) fits have rank(A) = #unknowns (unique solutions: 606 and 79 unknowns).
- **NUMERICAL** (`step_10.py`): B4 evaluated from its own HPL Taylor series to order 7000 — probing
  coefficients far beyond the fitted c^1280 — vs `corner_series` M4 at six angles (≈90-digit reference):
  |diff| ≤ 4.0·10⁻⁹² (π/6), 2.6·10⁻⁹² (π/4), 2.9·10⁻⁹² (π/3), 1.8·10⁻⁹³ (π/2), 4.9·10⁻⁹³ (2π/3),
  1.6·10⁻⁹³ (1.3 rad). At c = 1/√2 vs the 800-digit M4: **9.6·10⁻¹²¹**.

## 3. Task 3: depth structure — TO BE FILLED

## 4. Task 4: the hexagon corner closed form — TO BE FILLED
'''
with open(cwd + "/result.md", "w") as f:
    f.write(md)
print("result.md written, chars:", len(md))

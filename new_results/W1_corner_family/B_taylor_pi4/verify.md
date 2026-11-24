# W1-B verification: flat-angle Taylor bookkeeping for the π⁴ constant of Conjecture 4

**Verdict: CONFIRMED WITH CORRECTIONS.**
- The solver's negative result holds. The obstruction to Route B is real and exact, and I reproduced it with my
  own code. Every PROVED by-product I tested holds.
- The solver does **not** claim c₀ = c₁ = 0, and nothing here proves it. Conjecture 4 stays PROVED-MODULO-π⁴, with
  |c₀| < 10⁻³²² and |c₁| < 10⁻³²¹ certified (A3, A3b_certified2).
- The corrections are about scope and wording (§4). No number is wrong.

Working folder: `<scratch>/verify_W1_B_taylor_pi4/`. It holds a copy of the solver's folder, which I
re-ran, plus my own `v_*.py` scripts. Durable copies are in `B_taylor_pi4/verify/`: my scripts and their `.out`
files, and `verify/rerun/` with the re-run outputs of the solver's scripts. Every number below was printed by one
of these scripts. Labels: **PROVED**, **NUMERICAL**, **CONJECTURE**, as in the project.

---

## 1. Re-run of the solver's scripts

I copied the folder to scratch and ran all 12 scripts in the documented order: `hpl`, `zsum`, `mu3_symbolic 320`,
`b3_allN_proof`, `proof_controls`, `w4_taylor_formula`, `b2b2_test`, `obstruction_structure`,
`function_symbol_check`, `corner_equals_squircle`, `natural_extension_test`, `check_b3_split`.
- **Every archived `.out` file is reproduced byte for byte.** The two self-test scripts, which have no archived
  output, print `zsum self-tests passed: True` and HPL sanity residuals ≤ 2.3e-41.
- I found no hidden numerical step in the scripts labelled PROVED. Their identities are exact (sympy or Fractions),
  and none of them uses a float or PSLQ.
- There is no circularity: no step assumes Conjecture 4 or uses R. The only place R enters is
  `natural_extension_test.py`, which is a NUMERICAL illustration and is labelled as one.

## 2. Independent checks (my own code; none of it imports the solver's modules)

### 2.1 μ_{k,n} for symbolic n (`v_mu3_derive.py`): CONFIRMED, PROVED

- **Method.** This route differs from the solver's.
  - The u-expansion is an explicit binomial series (1+yW)^{1−2n} = Σ_j binom(1−2n, j) y^j W^j, with
    W = ((1−y)/y)^u − 1, times a truncated exp(uA). There is no series of a logarithm.
  - The Beta-derivative ratios come from sympy differentiating exp(logΓ(p)+logΓ(q)−logΓ(p+q)). There is no
    hand-written cumulant recursion.
  - Polygammas at n+j and 2n+j are rewritten with the shift rule and the values at integers.
- **The y-form of J_n.** I re-derived it by hand, including the g ↦ 1/g invariance (both exponents equal
  u−1+n(1−u)) and the y = g/(1+g) substitution. The prefactor C(2n,n)B(n,n)/2 = 1/n is also by hand.
- **Results:**
  - Euler's γ cancels for every k.
  - My μ_{3,n} minus the 22-monomial formula typed from result.md §2 is **0 identically in n**.
  - The same holds for G1's μ_{2,n} formula.
  - For k = 0..3 my formula equals G1's exact table `flat_exact_N320_K3.json` **for every n ≤ 320**.
- **Remark.** The solver's "independent exact evaluation" (n ≤ 40) shares its u-series and moment recursion with
  the symbolic code, so it is only a consistency check. My derivation is the independent check.

### 2.2 B3 and B2 for all n (`v_nested_proof.py`): CONFIRMED, PROVED

- **Method.** I wrote my own exact nested-sum algebra:
  - strict alternating Z-sums at 2n;
  - the quasi-shuffle product;
  - peeling shifts Z(2n+d) → Z(2n);
  - the HPL Taylor-coefficient recursion derived from f_a = 1/(1−at);
  - H^{(r)}_{n−1} = 2^{r−1}(Z_r+Z_{−r})(2n) − n^{−r}.

  The coefficients are my own exact Q(n) class built on python-flint `fmpq_poly`, not sympy. The μ side is taken
  from my own derivation (§2.1). G1's fit JSONs are parsed by my own parser.
- **Self-tests:** every rule is checked exactly at n = 1, 2, 3, 5, 8 (HPL coefficients also at n = 9 and 15).
  Result: True.
- **Results:**

  | Form | Part | Terms | Canonical words (B side / μ side) | [c^{2n}] B − μ ≡ 0 | [c^{2n+1}] B ≡ 0 | c⁻²…c³ exact |
  |---|---|---|---|---|---|---|
  | B3 | rational | 74 | 12/12 | True | True | True |
  | B3 | ζ(2) | 15 | 2/2 | True | True | True |
  | B3 | ζ(3) | 5 | 1/1 | True | True | True |
  | B2 | rational | 24 | 5/5 | True | True | True |
  | B2 | ζ(2) | 5 | 1/1 | True | True | True |

  μ_{2,n} has no ζ(3) part. The canonical forms also match direct values at n = 1, 2, 4, 7.
- **Controls:** I perturbed a random coefficient by +1/7, five times in R₃ and five times in Ψ. The identity
  failed in 5 of 5 trials each time.
- **Validity for small n:**
  - Every denominator is 2n+d with −1 ≤ d ≤ 3, n, n+1 or 2n±1, so none vanishes for n ≥ 1.
  - The peeling identity holds for every M ≥ 1. At n = 1 it reaches Z_v(0), which is correct.
  - The constant and c-linear terms sit only at c⁰ and c¹, which the exact low-order check covers.
- **Conclusion.** M₃ = B3 and M₂ = B2 (HPL forms) are **theorems** on |c| < 1. The only inherited caveat is
  G1-verify §3.1(ii), the one-sentence uniformity argument in Theorem A.

### 2.3 Taylor coefficients of W₄(c) (`v_w4_formula.py`): CONFIRMED, PROVED

- I mapped the solver's single-sum formula for w_n (typed from §4, with A = Z₋₁, B = Z₋₂, C = Z₋₃, P = Z₂, Q = Z₃)
  into Z-words with my quasi-shuffle. It equals my canonical form of [c^{2n}]W₄ **identically in n**, for all three
  parts. The solver checked this only for n ≤ 60, so it is now symbolic.
- The odd coefficients vanish identically.
- Directly from the HPL series:
  - w₀ = ½ζ(2) + ½ζ(3);
  - w₁ = −13/72 − (13/36)ζ(2) − ζ(3)/3;
  - w₂ = 2357/51840 − (109/1800)ζ(2) − ζ(3)/15.

  All three as in result.md.

### 2.4 The classicality obstruction (`v_symbol_obstruction.py`): CONFIRMED, PROVED

- **Symbols.** I compute the function-level symbol directly from dH_{a,w} = f_a H_w:
  S(H_{a₁…a_k}) = ℓ_{a_k} ⊗ … ⊗ ℓ_{a₁}, with ℓ₀ = [c], ℓ₁ = −[1−c], ℓ₋₁ = [1+c]. This is not Goncharov's formula,
  which the solver uses. The operator is 𝒪(a⊗b⊗c⊗d) = (a∧b)∧(c∧d) ∈ Λ²(Λ²V).
- **Controls** (200 random trials each). 𝒪 annihilates:
  - every Li₄-type symbol g⊗f⊗f⊗f;
  - every (1,3), (2,2) and (1,1,2) shuffle product.

  Result: True. This is the only property the theorem needs.
- **Words.** At weight 4 only the prefactors c⁻¹, c⁰, c¹ occur, with 16, 16 and 24 words (24 distinct).
- **Results:**
  - 𝒪(A₋₁) = −e₁₂∧e₂₃ − e₁₃∧e₂₃ and 𝒪(A₀) = e₁₂∧e₂₃ − e₁₃∧e₂₃, where e₁₂ = [c]∧[1−c], e₁₃ = [c]∧[1+c] and
    e₂₃ = [1−c]∧[1+c]. 𝒪(A₁) = 0.
  - Both nonzero forms equal the solver's ρ-forms exactly: −½ρ{−(1−c)/(1+c)}∧ρ{c²} and
    ρ{−(1−c)/(1+c)}∧(ρ{c}−ρ{−c}).
  - 𝒪(A₋₁) and 𝒪(A₀) are linearly independent (rank 2).
- **Reflection.** On the τ alphabet {τ, 1−τ, 1+τ, 1+τ²}, with constants dropped, the obstructions of A₋₁(c),
  A₀(c), A₋₁(s) and A₀(s) have **rank 4**.
  - The factorisations of c, 1±c, s and 1±s in τ are exact.
  - The map ι: τ ↦ (1−τ)/(1+τ) swaps c and s and fixes τ = √2−1 (sympy).
- **Specialisation at c₀ = 1/√2.** The letters 2, c, 1−c, 1+c lie exactly in ⟨√2, 1+√2⟩. The images of ρ{c},
  ρ{−c}, ρ{c²}, ρ{(1−c)/2}, ρ{−(1−c)/(1+c)} and ρ{(1−c)/(1+c)} are 1, −1, 0, −6, 6 and 4 times e = √2∧a.
  This matches the table in §6.
- **The theorem's logic** (result.md §5). I checked it step by step.
  - Steps 1, 3 and 4 are now exact computations.
  - Step 2, linear independence over the coefficient field, holds on any finite cover. The
    Deneufchâtel–Duchamp–Minh–Solomon criterion needs only that no nonzero ℂ-combination of the dlog forms is
    exact in the algebraic function field K. That follows from residues. Pick the letters multiplicatively
    independent modulo constants. An exact form has zero residues, and that forces Σa_j div(x_j) = 0, hence every
    a_j = 0.
  - Expanding the prefactors in a ℂ-basis of K that contains c⁻¹, 1 and c then isolates each A_e.
  - Algebraic prefactors are covered by the same argument. Rank 2 (rank 4 with the reflection) over ℚ is also
    rank 2 (rank 4) over ℂ, so no combination with algebraic-function prefactors is classical.
  - (C(c)^×/C^×)⊗ℚ → (K^×/C^×)⊗ℚ is injective, so 𝒪 stays nonzero after any algebraic extension of the arguments.

### 2.5 The square corner (`v_numeric_checks.py`, my own mpmath HPL code)

- **Lemma 1 (PROVED, by hand).**
  - With t = g/(1+g) the squircle integrand becomes the one stated in result.md §1.
  - The difference is ∫g^{u−1}(1+g)^{−1−u}dg + ∫(1+g)^{−1−u}dg − ∫g^{u−1}dg = 2^{−u}/u + (1−2^{−u})/u − 1/u = 0
    for u > 0, and for all u by analytic continuation. My overall sign is the opposite of the solver's, which does
    not matter since the result is 0.
  - With the note's eq. (Mk), M = −F/2, this gives M_k(π/2) = M_k.
- **W₂ (PROVED).**
  - The weight-2 HPL part of R₃ is exactly −H_{0,1} + H_{0,−1} = −½Li₂(c²), and Ψ contributes ½ζ(2).
  - So W₂(c) = ½ζ(2) − ½Li₂(c²) as functions, and W₂(c₀) = T exactly. The printed difference
    ½ζ(2) − ½Li₂(½) − T is 0.0.
- **Comparisons at c₀ (NUMERICAL, 185-digit working precision):**

  | Comparison | difference |
  |---|---|
  | B3(c₀) − M₃ (method A) | 8.6e-188 |
  | B3(c₀) − M₃ (method B) | 8.6e-188 |
  | W₂(c₀) − T | −9.2e-187 |
  | W₃′(c₀) − W₃′_note | −3.7e-186 |
  | W₄(c₀) − W₄_note | 6.4e-186 |

  The solver had 1e-130 for these. W₃′_note was re-read from main.tex Theorem 3.
- W₄(c₀) = 0.795785059832240811320359819629726762546527413871021247250728…, which agrees with main.tex.
- **Notation.** The values result.md quotes for A₋₁(c₀), A₀(c₀) and A₁(c₀) (1.0593, −1.0164, 0.4443) include the
  ζ(2) and ζ(3) pieces with the same prefactor. The pure weight-4 HPL parts are −0.16859, −1.01643 and 1.67217.
  The symbol statements concern the pure parts; the ζ-pieces have length ≤ 2 and play no role there.

## 3. Status of each claim

| Claim (result.md) | Status after verification |
|---|---|
| Route B obstructed: no classical R(c), with algebraic arguments and algebraic coefficients plus weight ≤ 3 terms, equals W₄(c) or M₃(α) on an interval | **PROVED** (my independent symbol computation; standard independence theorem) |
| δ₂,₂(A₋₁), δ₂,₂(A₀) ≠ 0 and δ₂,₂(A₁) = 0 (function level), with the explicit ρ-forms | **PROVED** (exact; reproduced) |
| Reflection α ↦ π−α: rank 4 | **PROVED** (exact; reproduced) |
| At c₀ all letters lie in ⟨√2, 1+√2⟩, so Λ²(Λ²) = 0; ρ table | **PROVED** (exact; reproduced) |
| (1) closed form of μ₃,ₙ for all n | **PROVED** (independent derivation; G1 table for n ≤ 320) |
| (2) B3 and B2 (HPL forms) hold for all n | **PROVED** (independent nested-sum proof) |
| (3) w_n closed form; no ζ(4); w₀ = ½ζ(2) + ½ζ(3) | **PROVED** (now symbolic in n) |
| (4) M(u;π/2) = squircle M(u) | **PROVED** (by hand) |
| (5) W₄ = B3(c₀) − W₃′ − T exactly; c₀ + √2c₁ = [π⁴-coefficient of an exact reduction] − (103/384 − 4963√2/5760) | **PROVED**, but a reformulation (see §4) |
| W₃′(c₀) = W₃′_note | **NUMERICAL** (now 3.7e-186) |
| Natural-extension table (§7) | **NUMERICAL** (re-run identical) |
| Two-variable lift Θ(y,z) (§8.4) | **CONJECTURE**, untested, correctly labelled |

## 4. Corrections

1. **Scope of the headline.** "Flat-angle Taylor bookkeeping cannot prove c₀ = c₁ = 0" should read: *no classical
   function-level extension of R exists along the corner family, with algebraic arguments and algebraic
   prefactors, nor along any reshuffling with the reflected family; so the Taylor comparison of step (3) has no
   valid input.*
   - That is what is proved.
   - It does not exclude other uses of the flat-angle data, for example Taylor control of the unobstructed
     component cA₁(c) combined with a point-level evaluation of A₋₁(c₀) and A₀(c₀).
2. **"B2 is a theorem" means the HPL form B2.**
   - G1's uniform Li₃ form B2′ (Li₃ at 8 arguments, ζ(3) coefficient 7/4 + 2c) still rests on G1's 10 HPL→Li₃
     dictionary identities. Those were checked only numerically (G1 `classical_dict.out`).
   - G1's open question 1 is settled for the HPL forms B1–B3, not for the classical rewriting B2′.
   - The ledger entry "Theorem B2 … with Li₃ at 8 standard arguments" should keep this distinction.
3. **By-product (5) carries no new information about c₀, c₁.**
   - It restates W₄ − R = (c₀ + √2c₁)π⁴ with W₄ written in HPL values at 1/√2. "No explicit π⁴" is true but
     uninformative: π⁴ is hidden in the weight-4 values H_w(1/√2).
   - The refinement that only W₄(c₀) remains uses the NUMERICAL W₃′ identity.
   - Uniqueness of "the" π⁴ coefficient of a reduction would need the unproved linear independence of R's basis.
4. **δ₂,₂ normalisation.** The displayed equalities "δ₂,₂(A_e) = …" hold in Λ²(Λ²) under the unnormalised map
   (a∧b)∧(c∧d), that is, after δ∧δ. Only non-vanishing is used, so nothing changes.
5. **Inherited caveats** (G1 verify). The split of μ₃,ₙ into rational, ζ(2) and ζ(3) parts is formal (the
   ℚ-independence of 1, ζ(2), ζ(3) is open). Theorem A needs one sentence on uniform convergence on a complex
   u-disc. Neither affects the results here.

## 5. Bottom line for Conjecture 4

Route B fails for a precise, exact reason, and I confirm it. The weight-4 part of M₃(α) has the non-classical
classes ρ{−tan²(α/4)}∧ρ{cos²(α/2)} (for the c⁻¹ prefactor) and ρ{−tan²(α/4)}∧(ρ{c}−ρ{−c}) (for c⁰). Both die only
at c₀, where the S-unit rank is 2.

c₀, c₁ remain rational numbers that are certified to be tiny but are not proved to be zero. Closing them needs a
point-level method:
- exact Kummer-type Li₄ relations at S-units of ℤ[√2, ½];
- a level-8 coloured-MZV reduction;
- or a denominator bound.

## 6. Files (`B_taylor_pi4/verify/`)

| Script | Checks | Label |
|---|---|---|
| `v_mu3_derive.py` (+`.out`, `v_mu_formula.pkl`) | independent symbolic μ_{k,n}; vs the typed formula and the G1 table (n ≤ 320) | PROVED |
| `v_nested_proof.py` (+`.out`) | own Z-sum/Q(n) engine; B3 and B2 for all n; mutation controls; W₄ canonical form | PROVED |
| `v_w4_formula.py` (+`.out`) | w_n single-sum formula vs canonical form, symbolic in n; w₀, w₁, w₂ | PROVED |
| `v_symbol_obstruction.py` (+`.out`) | own function symbols; 𝒪(A_e); ρ-forms; τ-rank 4; c₀ specialisation; controls | PROVED (exact) |
| `v_numeric_checks.py` (+`.out`) | exact W₂; B3(c₀) and its weight split vs the note at 185 digits | PROVED (W₂) / NUMERICAL |
| `rerun/*.out` | re-run outputs of all 12 solver scripts (all identical to the archived ones) | — |

To run: from `verify/`, run `v_mu3_derive.py 320` first, because it writes the pkl the others need. Then run the
other scripts in any order. Each takes seconds, except `v_mu3_derive.py 320`, which takes a few minutes.

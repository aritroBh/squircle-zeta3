# W1-B: flat-angle Taylor bookkeeping for the π⁴ constant of Conjecture 4

Folder: `research/10_frontier/W1_corner_family/B_taylor_pi4/`. Python: `.venv/bin/python`. Every number below is
printed by a script in this folder; each `.out` file is named after its script. G1's fitted HPL coefficients
(`G1_zeta3_any_angle/hpl_fit_k{2,3}_*.json`) and exact table (`flat_exact_N320_K3.json`) were read as data only.
Nothing outside this folder was modified.

Labels: **PROVED** means exact algebra or a complete argument given here. **NUMERICAL** means high-precision
agreement, with the digit count given. **CONJECTURE** covers everything else.

Notation: c = cos(α/2), s = sin(α/2), τ = tan(α/4), c₀ = 1/√2 (the square corner), b = √2−1, a = 1+√2.
The corner function is M(u;α) = Σ_k M_k(α)u^k, as in G1. HPLs H_w(c) use letters {0, 1, −1}
(Remiddi–Vermaseren).

---

## 0. Verdict

**NEGATIVE, with a precise obstruction (PROVED).** Route B cannot prove c₀ = c₁ = 0.

- **Obstruction.** The weight-4 part of the corner coefficient M₃(α), viewed as a function of c, is **not a
  classical polylogarithm**.
  - Its c⁻¹ and c⁰ components have a nonzero Λ²B₂ component of the cobracket.
  - So no extension R(c) of the right-hand side of conj:W4 can equal it on any interval. This holds whatever
    arguments (algebraic functions of c) and coefficient functions are chosen, and whatever weight ≤ 3
    corrections are added.
  - Step (3) of the route (comparing Taylor coefficients of R(c) with those of W₄(c)) therefore has no valid
    input.
- **Why the identity can still hold at c₀ = 1/√2.** At that point every symbol letter falls into the rank-2
  S-unit group ⟨√2, 1+√2⟩. There Λ²(O_S^×) is a line, so Λ²(Λ²) = 0, and the obstruction vanishes identically.
  - So classicality of W₄ is a property of the point c₀ alone. The flat-angle expansion cannot see it.
  - This matches the A3 Codex finding (reported in `A3_W4_proof/codex_b.log`, not re-run here) that the other
    natural deformation, the s-family, is also non-classical.

**By-products (all PROVED):**
1. The closed form of μ₃,ₙ for every n: rational, ζ(2) and ζ(3) parts (§2).
2. G1's HPL closed forms **B3 (M₃) and B2 (M₂) hold for all n**, hence as identities of analytic functions on
   |c| < 1 (§3). Before this they were EXACT through c⁶⁴⁰ (B3) and c⁶⁰⁰ (B2), plus NUMERICAL.
3. The closed form of the flat-angle Taylor coefficients of the weight-4 component W₄(c) (§4).
4. M(u; π/2) equals the squircle's M(u) exactly (§1). So the corner M_k at π/2 *is* the squircle M_k.
5. **π⁴ bookkeeping.**
   - W₄(c) contains no ζ(4) at the function level. Its Taylor coefficients involve only ℚ, ζ(2) and ζ(3).
   - Hence the proved representation W₄ = B3(c₀) − W₃′ − T contains no explicit π⁴ (and no Li₄(½)).
   - So c₀ + √2c₁ is the π⁴ produced by an exact reduction of the HPL values H_w(1/√2) to the basis of R, minus
     R's explicit π⁴ coefficient 103/384 − 4963√2/5760 (§8). Only 24 distinct weight-4 words occur.

---

## 1. Corner = squircle at α = π/2 (PROVED)

**Lemma 1.** M(u; π/2) = −F(u)/2, where F is the squircle's function (note, eq. (Mk)). Hence M_k(π/2) = M_k for
all k.

*Proof.* F = 2∫₀^{1/2}Φ dt. Put t = g/(1+g). Then
F = 2∫₀¹(1+g)^{−1−u}[g^{u−1} + 1 − g^{u−1}√(1+g^{2−2u})] dg.

Subtracting the corner integrand leaves
∫₀¹ g^{u−1}[1 − (1+g)^{−1−u}] dg − ∫₀¹(1+g)^{−1−u} dg = (1 − 2^{−u})/u − (1 − 2^{−u})/u = 0.
This uses ∫₀¹g^{u−1}(1+g)^{−1−u}dg = 2^{−u}/u. ∎

Check (`corner_equals_squircle.out`): the two integrals agree at u = 0.3, 0.7 and −0.4, to 40 digits
(1e-27 at u = −0.4, where the integrand is near-singular).

## 2. Task (1): exact μ₃,ₙ for symbolic n (PROVED)

Script: `mu3_symbolic.py`, output `mu3_symbolic.out`. This is G1's Theorem A re-implemented independently, as
sparse polynomial algebra over ℚ(n). It uses the moments D_{i,l} of the Beta-function cumulants, the polygamma
values at p = n+j, q = n, and the overall factor −1/(n(2n−1)).

Write D_n = H_{2n−1} − H_{n−1}, h_k = H^{(k)}_{n−1} and g_k = H^{(k)}_{2n−1}. Then for every n ≥ 1:

```
mu_{3,n} = - z3/(4n^2-1)
           + z2 [ (5n^2+5n+1)/(2n(n+1)(2n+1)^2) - D_n/(4n^2-1) ]
           + D_n^3/(6n(2n-1)) - D_n^2/(4n^2) + D_n g2/(2n(2n-1)) - D_n h2/(2n(4n^2-1))
           + (4n^2-2n-1) D_n/(4n^3(2n+1)) + g3/(3n(2n-1)) + (n-1) h3/(3n(4n^2-1))
           - g2/(4n^2) - (3n+1)(2n^2-1) h2/(4n^2(n+1)(2n+1)^2)
           - (16n^5-12n^3+4n^2+5n+1)/(8n^4(n+1)(2n+1)^2)
```

- The script prints these 22 monomials. The regrouping into powers of D_n is read off the printed coefficients;
  for example, the h1³, h1²g1, h1g1² and g1³ terms combine to D_n³/(6n(2n−1)).
- Example: μ₃,₁ = 35/72 − ζ(2)/36 − ζ(3)/3.
- Checks:
  - The formula equals an independent exact evaluation (numeric polygamma values, Fractions) for n ≤ 40.
  - It equals G1's stored exact table for all n ≤ 320.
  - The same holds for k = 0, 1, 2. The k = 2 output reproduces G1's μ₂,ₙ.
- **New relative to G1:** the closed forms of the ζ(2) part and the rational part of μ₃,ₙ for symbolic n. G1 had
  the ζ(3) coefficient −1/(4n²−1) symbolically and the closed form of μ₂,ₙ. The other parts of μ₃,ₙ existed only
  as exact tables (n ≤ 320) and through the fitted B3.
- Caveats, carried over from G1's verifier:
  - (i) The split into rational, ζ(2) and ζ(3) parts is the canonical *formal* split from the polygamma
    reduction. Uniqueness for the real number μ₃,ₙ would need the open ℚ-independence of 1, ζ(2), ζ(3). Nothing
    below depends on that uniqueness.
  - (ii) Theorem A needs the one-line uniformity argument for termwise u-differentiation (verify.md §3.1).

## 3. B3 and B2 for all n (PROVED)

Script: `b3_allN_proof.py`, output `b3_allN_proof.out`. Engine: `zsum.py`; self-tests are in its `__main__`.

**Method.**
- Both sides are rewritten exactly as ℚ(n)-linear combinations of alternating Z-sums Z_w(2n).
- **μ side:** H^{(r)}_{n−1} = 2^{r−1}(Z_r + Z_{−r})(2n) − n^{−r} and H^{(r)}_{2n−1} = Z_r(2n) − (2n)^{−r}.
  Products are expanded with the quasi-shuffle.
- **HPL side:** [c^N]c^eH_w = h_w(N−e), with h_w(N) = Σ coef·ε^N N^{−k} Z_v(N−1) computed exactly from the HPL
  recursion.
- At N = 2n or 2n+1 the signs ε^N are constants. The shifts Z_v(2n+j), with −2 ≤ j ≤ 2, reduce exactly to
  Z-sums at 2n. All identities used are valid for n ≥ 1.
- Equal canonical forms are sufficient for equality; no independence statement is needed.

**Results.** For each of the three parts of B3 (rational R₃, 74 HPL terms; ζ(2) part Ψ, 15 terms; ζ(3) part
−Φ, 5 terms):
- `[c^{2n}]B3 − μ₃,ₙ ≡ 0` identically in n;
- `[c^{2n+1}]B3 ≡ 0` identically in n;
- the coefficients c⁻², …, c³ agree exactly.

The script prints **"B3 PROVED for all n (all three parts): True"**. The same check run on G1's B2 (the M₂ form:
24 + 5 terms) prints "B2 (HPL form of M_2) PROVED for all n: True" (`proof_controls.out`).

**Controls** (`proof_controls.out`):
- Changing one B3 coefficient by +1/6 makes both identities fail.
- The symbol code is checked three ways: shuffle compatibility on 136 word pairs, the x → −x symmetry for all
  words of weight ≤ 4, and the last-slot derivative identity for all words of weight 2–4.

So G1's open question 1 is settled for M₂ and M₃: **B2 and B3 are theorems.**

## 4. The weight-4 component W₄(c) (formula PROVED; match with the squircle W₄ NUMERICAL)

**Definition.** Take the weight-4 terms of B3:
W₄(c) = [weight-4 HPL terms of R₃] + ζ(2)·[weight-2 part of Ψ] + ζ(3)·[weight-1 part of −Φ]
= c⁻¹A₋₁(c) + A₀(c) + c·A₁(c).

The pure parts A_e use 16, 16 and 24 HPL words of weight 4 (24 distinct words in all; `corner_equals_squircle.out`).

**Taylor coefficients (PROVED).** W₄(c) = Σ_{n≥0} w_n c^{2n}, with w₀ = ½ζ(2) + ½ζ(3). For n ≥ 1, with
A = H_n − H_{2n}, B = ½H_n^{(2)} − H_{2n}^{(2)}, C = ¼H_n^{(3)} − H_{2n}^{(3)}, P = H_{2n}^{(2)} and
Q = H_{2n}^{(3)}:

```
w_n = - z3/(4n^2-1) + z2 [ A/(4n^2-1) - (6n+1)/(4n(2n-1)(2n+1)^2) ]
      - (7n-4)/(4n^4(2n-1)) + (3n-2)A/(2n^3(2n-1)) + (8n^3+6n^2-5n-2)B/(2n^2(2n-1)(2n+1)^2)
      + (6n+1)P/(2n(2n-1)(2n+1)^2) + 2Q/(4n^2-1) - C/(n(4n^2-1)) + (AB-Q)/(n(4n^2-1))
      - (n-1)(A^2-P)/(2n^2(2n-1)) - 2(AP-C)/(4n^2-1) - (A^3-3AP+2C)/(6n(2n-1))
```

- Where this comes from: `b3_allN_proof.out` gives the canonical Z-sum form. `w4_taylor_formula.out` checks
  that it equals the single-sum form above (n ≤ 60) and the exact Taylor coefficients of the weight-4 part of B3
  (n ≤ 150). It also checks that the odd coefficients vanish through c²⁹⁹.
- Examples: w₁ = −13/72 − (13/36)ζ(2) − ζ(3)/3 and w₂ = 2357/51840 − (109/1800)ζ(2) − ζ(3)/15.
- The ζ(3) part of w_n equals the whole ζ(3) part of μ₃,ₙ. **No ζ(4) appears in any w_n** (the "π⁴
  bookkeeping" of the brief).

**Match with the squircle (NUMERICAL).** By Lemma 1, B3 at c₀ equals the squircle M₃ exactly. Numerically
(`check_b3_split.out`, 130-digit working precision), the weight split of B3 at c₀ reproduces the squircle split
term by term:

| Comparison at c₀ | difference |
|---|---|
| B3(c₀) − M₃ (1200-digit method A) | 2.8e-130 |
| W₂(c₀) − T | 1.6e-130 |
| W₃′(c₀) − W₃′_note | 3.6e-130 |
| W₄(c₀) − W₄_note | −2.7e-130 |

Values: W₄(c₀) = 0.795785059832240811320359819629726762546527413871021247250728…, with
A₋₁(c₀) = 1.0592837626503195519…, A₀(c₀) = −1.0164301034434735106… and A₁(c₀) = 0.44429173654204894027…

The exact identity W₂(c₀) = T is proved (½ζ(2) − ½Li₂(½) = (π²+6λ²)/24). The weight-3 identity W₃′(c₀) = W₃′_note
is NUMERICAL only. It is not needed below: the obstruction lives in the top-weight HPL coefficients of M₃
itself.

## 5. Tasks (2)–(3): no natural extension of R exists (PROVED)

**Theorem (obstruction).** Let I ⊂ (0,1) be an interval. There is no identity on I of the form

  M₃(α) = F(c) + G(c),

where:
- F is a finite sum Σ p_i(c)Φ_i(c), with p_i algebraic functions of c and each Φ_i either Li₄(x(c)) or a
  product of classical polylogarithms Li_k(x(c)) and logarithms of total weight 4, with arguments x(c)
  algebraic over ℚ(c) and arbitrary complex coefficients;
- G is any combination of iterated integrals of weight ≤ 3 with algebraic coefficients, for example any
  constants times lower-weight polylogarithms.

In particular, no extension R(c) of the right-hand side of conj:W4 satisfies R(c) = W₄(c) near c₀. This holds for
any choice of the arguments (c/(1+c), (c−1)/c, −(1−c)/(1+c), c, 1−c, 1/(1+c), …) and of the coefficient functions
(√2 ↦ 1/c, 2c, …).

*Proof.*
1. By §3, M₃ = B3 exactly. Its weight-4 HPL part is Σ_e c^e A_e with rational coefficients. The ζ(2) and ζ(3)
   terms have length ≤ 2 as functions.
2. Iterated integrals of dlog forms are linearly independent over the coefficient function field
   (Deneufchâtel–Duchamp–Minh–Solomon 2011 for rational letters; the same argument applies on the finite cover
   where the algebraic arguments live). So an identity of the stated form would force each A_e to be a
   ℂ-combination of length-4 parts of classical pieces, after expanding the p_i in a basis containing c⁻¹, 1, c.
3. Every such piece has vanishing Λ²B₂ cobracket δ₂,₂ (Goncharov; Goncharov–Spradlin–Vergu–Volovich 2010):
   - δ{x}₄ = {x}₃⊗x;
   - weight-1 × weight-3 products die under the inner wedge;
   - products of two weight-2 factors die under the outer wedge.
   δ₂,₂ is natural under extension of the function field (ℚ-vector spaces of letters inject).
4. The exact computation gives δ₂,₂(A₋₁) ≠ 0 and δ₂,₂(A₀) ≠ 0. ∎

**The obstruction, explicitly** (`b2b2_test.out`, `obstruction_structure.out`). Write ρ{z} := z∧(1−z) for the
Λ²-image of {z}₂. Then:

  δ₂,₂(A₋₁) = −½ ρ{−(1−c)/(1+c)} ∧ ρ{c²},
  δ₂,₂(A₀) = ρ{−(1−c)/(1+c)} ∧ (ρ{c} − ρ{−c}),
  δ₂,₂(A₁) = 0.

So the c¹ component is not obstructed (δ₂,₂ = 0). Note that −(1−c)/(1+c) = −tan²(α/4). A single
dilogarithm class, {−tan²(α/4)}₂, carries the whole obstruction.

- **Function-level version** (`function_symbol_check.out`): with the constant letter 2 set to 0 (the symbol of
  W₄ as a function of c), δ₂,₂(A₋₁) and δ₂,₂(A₀) each still have 2 nonzero components, and δ₂,₂(A₁) is still 0.
  This is the version the theorem needs.
- **Reflection α ↦ π−α does not help.** In τ the map ι: τ ↦ (1−τ)/(1+τ) swaps c and s and fixes the square
  corner τ = b. In the τ-alphabet {τ, 1±τ, 1+τ²} the obstructions of A₋₁(c), A₀(c), A₋₁(s) and A₀(s) have
  **rank 4**, in both symbol versions.
- Hence no combination with rational-function prefactors is classical, including the reshufflings
  q(τ)[A_e(c) + A_e(s)] with q(b) = c₀^e that agree with W₄ at the fixed point.

## 6. Why classicality survives at c₀ (PROVED)

Specialise the letters at c₀: 2 = (√2)², c = (√2)⁻¹, 1−c = b/√2, 1+c = a/√2. Every class then becomes a multiple
of e = √2∧a, which spans the line Λ²(O_S^×)⊗ℚ:

| class | ρ{c} | ρ{−c} | ρ{c²} | ρ{(1−c)/2} | ρ{−(1−c)/(1+c)} | ρ{(1−c)/(1+c)} |
|---|---|---|---|---|---|---|
| image at c₀ | e | −e | 0 | −6e | 6e | 4e |

(`obstruction_structure.out`.)

- δ₂,₂(A₋₁) dies because ρ{c²} → ρ{½} = 0, that is, because Li₂(½) is a product.
- δ₂,₂(A₀) dies because 6e∧2e = 0.

This is the precise sense in which the identity of conj:W4 is pointwise: it uses the arithmetic of ℚ(√2) (rank-2
S-units, so Λ²(Λ²) = 0), and that arithmetic is invisible to any expansion in c.

## 7. Illustration: the literal natural extension (NUMERICAL)

`natural_extension_test.py` builds R_nat(c) by substituting:
- b → c/(1+c), −b → (c−1)/c, −b² → −(1−c)/(1+c);
- 1/√2 → c, 1−1/√2 → 1−c, 2−√2 → 1/(1+c);
- L → log((1+c)/c), and √2 → 1/c (variant A) or 2c (variant B).

Each substitution reduces correctly at c₀. Results at 60 digits:

| c | W₄(c) | R_nat − W₄ (A) | R_nat − W₄ (B) |
|---|---|---|---|
| 0.05 | 1.42055653223 | 1.38e5 | 1.86e3 |
| 0.5 | 1.12058887535 | 816.75 | 429.07 |
| 0.7 | 0.80929087723 | 15.66 | 15.33 |
| 1/√2 | 0.79578505983 | 1.3e-58 | 1.8e-58 |
| 0.72 | 0.77080697301 | −27.07 | −27.93 |

The literal extension agrees with W₄ only at c₀, where it *is* R, and it blows up toward the flat angle. The theorem
of §5 says that no other choice of arguments or coefficients can do better.

## 8. Consequences for c₀, c₁

1. **Flat-angle Taylor bookkeeping cannot determine c₀, c₁ (PROVED, §5).** The intended mechanism works only for
   classical pieces, where "no ζ(4) in the Taylor coefficients" would fix the π⁴ term. The obstructed components
   A₋₁ and A₀ cannot be put in that form along the corner family.
2. **What the bookkeeping does give (PROVED).**
   - Exactly, W₄ = B3(c₀) − W₃′ − T, by §3, Lemma 1 and the note's Theorem 3.
   - B3(c₀) is a ℚ(√2)-combination of HPL values H_w(1/√2) of weights 1–4 (24 distinct weight-4 words), with
     rational, ζ(2) and ζ(3) coefficients. W₃′ and T are the note's classical weight ≤ 3 constants.
   - **No explicit π⁴ (and no Li₄(½)) occurs anywhere.**
   - With the numerically verified weight-3 identity W₃′(c₀) = W₃′_note (§4), only the weight-4 component W₄(c₀)
     remains, and that component has no log 2 either.
   - So Conjecture 4 ⟺ this HPL expression equals R. Then c₀ + √2c₁ is the π⁴ that an exact, constant-tracking
     reduction of the H_w(1/√2) to R's Li₄ basis produces, minus 103/384 − 4963√2/5760.
   - The alphabet here is {0, ±1} plus the point 1/√2, which may be a more convenient starting point for an
     exact reduction than the six-point squircle alphabet.
3. **What is still needed.** A point-level exact method. Options:
   - exact Li₄ functional equations (Kummer type) at S-units of ℤ[√2, ½], with tracked constants;
   - an exact weight-4 reduction of level-8 coloured MZVs, since the Cayley map sends the squircle letters to
     8th roots of unity;
   - a denominator bound for c₀, c₁, combined with the certified |c_i| < 10⁻³²¹.

   Function-level coaction arguments along the corner family (this report) or the s-family (A3, codex_b.log)
   cannot work.
4. **Untested idea (CONJECTURE).** The obstruction is X(c)∧Y(c), where the class Y becomes ρ{½} = 0 (for A₋₁) or
   proportional to X (for A₀) at c₀.
   - A two-variable lift Θ(y,z) with Θ(c,c) = W₄(c) and δ₂,₂Θ = X(z)∧Y(y) would make z ↦ Θ(c₀,z) classical at
     the function level, so its constants could then be fixed by Taylor bookkeeping in z.
   - This needs the M₃ Euler-substitution integral with the letters in s and τ decoupled.
   - The last comparison with R would still need exact Li₄ relations, unless the arguments land in R's 7-element
     basis.
   - Not attempted here.

## 9. Scripts and outputs

| Script | What it does | Label |
|---|---|---|
| `hpl.py` | HPL Taylor coefficients (exact and mpmath), values, constant-retaining symbols | tool (self-tests) |
| `zsum.py` | alternating Z-sum algebra over ℚ(n): quasi-shuffle, shifts, harmonic numbers, HPL coefficients | tool (self-tests pass) |
| `mu3_symbolic.py` | closed form of μ_{k,n}, k ≤ 3; checks against independent evaluation and G1 table; Z-sum form | PROVED |
| `b3_allN_proof.py` | B3 for all n; Z-sum form of the W₄ Taylor coefficients | PROVED |
| `proof_controls.py` | mutation test; B2 for all n; symbol consistency checks | PROVED / control |
| `w4_taylor_formula.py` | single-sum form of the W₄ Taylor coefficients, exact checks | PROVED |
| `check_b3_split.py` | B3 and its weight split at c₀ vs the squircle values (130 digits) | NUMERICAL |
| `b2b2_test.py` | δ₂,₂ of A₋₁, A₀, A₁ (constant-retaining), product/Li₄ sanity checks | PROVED (exact) |
| `obstruction_structure.py` | δ₂,₂ in terms of ρ{z}; τ-pullback rank under α ↦ π−α; specialisation at c₀ | PROVED (exact) |
| `function_symbol_check.py` | the same with the constant letter dropped (function-level symbol) | PROVED (exact) |
| `natural_extension_test.py` | the literal extension R_nat(c) vs W₄(c) | NUMERICAL |
| `corner_equals_squircle.py` | numerical check of Lemma 1; HPL word count | NUMERICAL |

Run order: `mu3_symbolic.py` writes `mu_zsum.pkl`; then `b3_allN_proof.py` writes `w4_taylor_zsum.pkl`; then the
rest. Each script runs in seconds except `mu3_symbolic.py 320`, which takes a few minutes.

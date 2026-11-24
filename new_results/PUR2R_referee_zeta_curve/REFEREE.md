# PUR2R: referee report on PUR2 (the ζ-generating curve)

Under review: `../PUR2_zeta_generating_curve/result.md`. PUR1 and its REFEREE.md are taken as proved.
I did not reuse PUR2's scripts. Everything below comes from my own code in this folder.

**Verdict in one line.** All the mathematics is correct. Three statements need a wording repair (Lemma 1 at the
point α = β = 0; Theorem A at c = ±1; the support-function step, which cites Schneider "from memory"). One
column of a numerical table is a grid artefact. No claim is false.

| file | what it does |
|---|---|
| `r_curve.py` | my implementation of K_{τ,c}: Λ by safeguarded Newton (no Lambert W, no Cardano), h, h′, h+h″, closed form − tail |
| `r_numerics.py` → `r_numerics_output.txt` | perimeters at τ = 0.3, 0.9 by u-quadrature vs closed form − tail (dps 60); full-circle quadrature; Lemma 1 at complex α, β; c = 1 at τ = 1 |
| `r_convex.py` → `r_convex_output.txt` | h+h″ formula vs mp.diff; grid minima (4000 points); h′ jump; E1/E2 lower bounds; polygon reconstruction; spot checks of E2(iii) |
| `r_interval.py` → `r_interval_output.txt` | python-flint (arb, 200 bits) certificates at κ² = 9/7, and a bracket for κ_max |
| `r_symbolic.py` → `r_symbolic_output.txt` | (A) coefficient table and (C) by sympy; (D) by an independent route (moments γ_kb by quadrature, not the Beta reflection) |
| `r_extra.py` → `r_extra_output.txt` | Remark A′ for c = 2, past the pole τ = 2/3; κ² = 9/7 at τ = 1.5 |

---

## (A) Family formula: PROVED AFTER REPAIRS (wording only)

**Definition of the body.** Checked line by line.
- p(u) = 1 + sin u − cos u has p′ = M > 0 on [0, π/4], with p(0) = 0 and p(π/4) = 1.
- So Λ_c = Φ_c^{−1}∘p is well defined, smooth and strictly increasing exactly when Φ_c is a smooth increasing bijection [0, V] → [0, 1], that is, when φ_c > 0 on [0, V].
- Real c: φ_c > 0 everywhere, so Λ_c is always defined. Imaginary c: see (C).
- The body depends only on (τ, c), and only on c². c and −c give the same body.
- On each window, ∫ M f du = ∫₀^V φ (1+e^{−2v/τ})^τ dv and ∫₀^V φ = p(π/4) = 1. So ∫₀^{2π} h = 8 + 8F is correct.

**Lemma 1.** Correct, with one gap in step 2.
- The domain, the bound |(1+x)^τ − 1| ≤ |τ|2^{|τ|}x (via e^{τg} − 1 with g = log(1+x) ≤ x), and local uniform convergence are all right.
- Step 3 is right. On U both Re α < 0 and Re β < 0 hold. The split ½∫₀¹x^{−α−1}(1+x)^τ dx + 1/(2α), the inversion x ↦ 1/x, and the Beta parameters p = −α, q = −β are all correct.
- **Gap (step 2).** Step 2 removes the pole on α = 0 (β ≠ 0) and on β = 0 (α ≠ 0). It does not treat the point α = β = 0, which lies in D. There Γ(−τ) also has a pole.
- **Repaired text for step 2:** "Write Γ(−z) = −Γ(1−z)/z. The right side equals ½·(τ/(αβ))·(1 − R(α,β)), with R = Γ(1−α)Γ(1−β)/Γ(1−α−β), analytic on D. R(0,β) = R(α,0) = 1, so 1 − R = αβ·g(α,β) with g analytic on D. Hence the right side is ½τg, analytic on all of D."
- This repaired form also gives the "Rewriting" paragraph directly.
- Own check: α = 0.3+0.7i, β = −0.45−0.7i gives 2·10⁻⁴⁸. α = 0.6+0.9i, β = 0.2−0.4i gives 9·10⁻²⁸ (slow oscillatory decay limits the quadrature).

**Theorem A.** Correct.
- v = τs gives F_∞ = (τ/2)∫(e^{2aτs} + e^{2bτs})B ds, and Lemma 1 then gives F_∞ = (1−R)/(4ab).
- Tail: ∫_V^∞ e^{2av}e^{−2nv/τ}dv = τe^{(2a−2n/τ)V}/(2n−2aτ), so 8T = E_c as printed. Fubini holds because n/τ > Re a for all n ≥ 1.
- The range τ < 1/Re a is right: τ < 2/(1+|c|) for real c and τ < 2 for imaginary c.
- At c = 0 the formula reduces to PUR1 §6, and the leading tail term is right.
- **Repair (c = ±1).** Here ab = 0, so (2/ab)(1−R) is 0/0 as printed. **Add:** "For c = ±1 the bracket (2/(ab))(1−R) means its limit, which is analytic in c because 1 − R vanishes at ab = 0. Theorem B gives it explicitly."

**Remark A′ (continuation, pole cancellation).** PROVED.
- The left side is real-analytic in τ > 0.
- E_c converges locally uniformly in a complex neighbourhood of (0, ∞) away from n = aτ and n = bτ, because |e^{−2nV/τ}| decays geometrically and |C(τ,n)| grows only polynomially.
- So the right side is meromorphic there and equals an analytic function, and its poles must cancel.
- Own checks:
  - c = 2 (direct range τ < 2/3): agreement to 10⁻³⁰ at τ = 0.5, 0.8 and 1.
  - c = 2 at τ = 2/3 + 10⁻¹²: agreement to 2.6·10⁻¹⁹, as expected from cancellation near the pole at dps 30.
  - κ² = 9/7 at τ = 1.5: agreement to 0.

**Coefficients.**
- Newton's recursion, p_k − 1 for k ≤ 7, and the table of τ²…τ⁷ are all confirmed by sympy (`r_symbolic_output.txt`: "matches: True" for n = 4…7).
- At e = 0 the series is 2Σζ(n)τⁿ.

## (B) c = 1: PROVED

- e^{2Λ} + 2Λ = 1 + 4p, so E = e^{2Λ} satisfies E + log E = q and Λ = ½[q − W(e^q)]. Both are right.
- My Newton-computed V and Λ(0.37) agree with the Lambert form to 5·10⁻⁵⁷ and 4·10⁻⁶².
- Each step of the direct proof checks:
  - the reduction 8F_∞ = 2τI;
  - the inversion;
  - J(Y) via 1 + x = 1/(1−t), where (1+Y)^τ − Y^τ → 0 needs τ < 1;
  - the monotone convergence of Σ(t^{n−τ} − tⁿ);
  - ψ(1−τ) = −γ + Σ[1/(n+1) − 1/(n+1−τ)].
- The tail is Theorem A's with a = 1, b = 0.
- τ = 1 value: 8∫h − (8 + 4V + 2(1−e^{−2V})) = −2.6·10⁻⁵⁶. Closed form − tail at τ = 1 − ε differs from it by −0.753ε (ε = 10⁻¹⁰, 10⁻²⁰), so the formula is continuous at τ = 1.

## (C) c = iκ: PROVED

- Φ = [e^v(cos κv + κ sin κv) − 1]/(1+κ²) is correct.
- The admissibility criterion κe^{π/(2κ)} > 2 + κ² is correct. Strictness is needed: at κ_max, φ(V) = 0.
- The monotonicity proof is correct: q′ has the sign of −2κ³ − πκ² + 4κ − 2π < 0.
- **κ_max is now certified** (arb): q(1.50781081090595335) > 0 > q(1.50781081090595336).
- At κ² = 9/7 (arb): Φ(π/(2κ)) = 1.5448654… > 1 and V ∈ [0.77678030588986128137171322132206368973 ± 10⁻⁴⁰].
- The formula holds for every κ as an identity of functions. It is the perimeter of a curve only for κ < κ_max.
- The τ⁴ coefficient is (2 − 7e/2)ζ(4), and e = (1+κ²)/4 = 4/7 ⟺ κ² = 9/7. τ⁵ = (6/7)ζ(5) − (8/7)ζ(2)ζ(3) = −1.3709811816454499125…. τ⁶ = −(22/49)ζ(6) − (4/7)ζ(3)² = −1.2824467289755028494…. All three confirmed symbolically.

## (D) Universality: PROVED (and the proof can be made elementary)

I re-derived c₂ … c₇ by a route independent of the Beta reflection.
- PUR1 (a) gives F_n = Σ_{k+b=n−1} φ_k k! γ_kb, with φ_k the Taylor coefficients of e^vψ.
- I computed γ_kb by quadrature (dps 45).
- Coefficientwise in (1, ψ″, ψ⁗, ψ⁽⁶⁾), all of c₂ … c₇ as printed agree to ≤ 4·10⁻⁴⁴.

**A two-line proof of c₄ and c₅ that a paper may prefer.**
- F_n is linear in φ_k for k ≤ n − 2.
- With ψ even, φ_k for k ≤ 3 depends on ψ only through ψ″(0): φ₂ = ½ + ψ″/2 and φ₃ = 1/6 + ψ″/2. So c₄ and c₅ are affine in ψ″(0).
- Two values fix them: ψ″ = 0 (PUR1's square: (9/8)ζ(4) and (3/2)ζ(5) − ½ζ(2)ζ(3)) and ψ″ = 1 (Theorem B: 2ζ(4) and 2ζ(5)).
- This gives exactly c₄ = (9+7ψ″)ζ(4)/8 and c₅ = (3+ψ″)ζ(5)/2 + (ψ″−1)ζ(2)ζ(3)/2.
- The ψ″-coefficient of c₄ is 8γ₂₁ = η(4) = 7ζ(4)/8, as it must be.
- In the same way c₆ and c₇ are affine in (ψ″, ψ⁗), and the family cosh(cv) traces the non-degenerate curve (c², c⁴), so §2's table determines them. Their printed forms pass this test.

**Rigidity.**
- ψ^{(2k)}(0) enters c_{2k+2} only through φ_{2k}·(2k)!γ_{2k,1}, with coefficient 8γ_{2k,1} = 8η(2k+2)/2^{2k+1} = 4^{1−k}η(2k+2) > 0. Numerically this matches for k = 1, 2 to 30 digits.
- c_{2k+3} involves no new derivative.
- So the induction behind "c_n = 2ζ(n) for n ≤ N ⟺ ψ^{(2j)}(0) = 1 for 2j ≤ N − 2" is valid. The class assumes an even Taylor series, so ψ′(0) = 0 is already built in.
- The correction in PUR2 (ψ⁗(0) does not enter c₄ or c₅) is right.

## (E) Convexity

**Formula h + h″ = 2fσ(M/φ)B.** Re-derived: f′ = −2fσΛ′, σ′ = −σ(1−σ)·2Λ′/τ, Λ″ = M′/φ − M²φ′/φ³.
- My implementation agrees with mp.diff to relative error 3·10⁻³⁵.
- The split 2K = 2 + 2(1−σ)(1/τ−1) is right.
- The key identity checks: w = M′ = 1 − p (cos u − sin u = 1 − p), M² = 2 − w², dφ/dp = 1 + χ, and the constant part 6p.

**u = 0.** h′(0⁺) = f(0)[1 − 2·½·1] = 0. h(u) = H(|u|) with H smooth on [0, π/4], so h ∈ C² (not C³ in general). (h+h″)(0) = 2^τ(1/τ − 1). PROVED.

**u = π/4.** h′(π/4⁻) = √2·f′ = −4fσ/φ(V), and by h(π/2 − u) = h(u) the jump is 8fσ/φ(V) > 0.
- One-sided mp.diff agrees to 16 digits (limited by the one-sided difference). The analytic h′ agrees to 10⁻³⁵ at interior points.
- *Only two windows meet* at each diagonal π/4 + jπ/2, not four. h is continuous there, with m = √2 and Λ = V on both sides.

**Support-function step: PROVED AFTER REPAIR (citation).** The criterion they use is correct: a continuous 2π-periodic h is the support function of a planar convex body iff h + h″ ≥ 0 as a distribution. For piecewise-C² h, that holds iff it holds on each open piece and every jump of h′ is ≥ 0. Replace the citation "from memory" by this self-contained argument:
- Put x(u) = h cos u − h′ sin u and y(u) = h sin u + h′ cos u on each piece. Then (x′, y′) = (h+h″)(−sin u, cos u).
- So where h+h″ ≥ 0 the tangent direction is u + π/2, which is monotone.
- At a positive jump J of h′, the curve jumps by J·(−sin u, cos u), a segment in the same direction.
- The resulting closed curve has a monotone tangent angle with total turning 2π, so it bounds a convex body whose outer normal at x(u) is e_u. Hence h is its support function.
- Its length is ∫(h+h″) over the circle, atoms included. Since h′ is periodic, the integral of h″ (atoms included) is 0, so the length is ∫₀^{2π}h du.
- **So "perimeter = ∫h" is legitimate, short edges included.** My polygon test confirms it numerically.

**E1 (log ψ convex ⇒ convex for all τ ∈ (0,1]).** PROVED.
- ψ even gives χ(0) = 0, and log ψ convex makes χ ≥ 0 and nondecreasing, so ∫₀^pχ ≤ pχ(p).
- With w = 1 − p ≥ 0, (2 − w²) − 3wp = 1 − p + 2p², whose discriminant is 1 − 8 < 0. So Dφ ≥ 6p.
- log cosh(cv) is convex. Grid check for c = 1: min(Dφ − 6p) = 3.9·10⁻⁴ ≥ 0.

**E2(i), κ² = 9/7.** PROVED, re-certified independently.
- The chain checks: −3w∫χ ≥ 0, 2 − w² = 1 + 2p − p² ≤ 2, k′ = κ²sec²(κv)/φ, and the condition k′ ≤ 3 ⟺ κ² ≤ 3e^v cos³κv.
- The unimodality argument is right.
- **My arb certificate does not need unimodality.** Subdividing [0, V_hi] into 4000 balls gives min 3e^v cos³(κv) ≥ 1.68219… > 9/7, and at V the ball is 1.682519737300218665559563 ± 5·10⁻²⁶.
- PUR2's mpmath-iv certificate is also valid in logic: Φ is monotone, so the bracket on Φ(V ∓ 10⁻²⁰) encloses V. arb is the more trustworthy library.

**E2(ii).** PROVED.
- 2(1−σ) ≥ 1 because σ ≤ ½ (this needs τ ≤ 1).
- The term 3wI/(2−w²) ≤ 0 because I ≤ 0, and −χ ≤ κ tan κV.
- The general version (any ψ, τ ≤ 1/(1+δ)) is also right.

**E2(iii).** NUMERICAL, as labelled. Spot checks: κ = 1.3 gives τ* = 1, and κ = 1.4 gives τ₀ = 0.2084, τ* = 0.6982 (grid 600), matching the table.
- **Grid artefact.** For κ ≤ κ_D, and in the real-c list of §6, the "min Dφ/M²" column is not a true minimum. D(0) = 0, so the infimum is 0, approached as u → 0. The printed positive numbers are just the value at the first grid point: I get 0.0056 at κ = 1.3 on a grid of 600 against their 0.011 on 300.
- Only the sign of D near π/4 is meaningful there. Reword the column or drop it.

**E3 (τ > 1 never convex).** PROVED. D(0) = 0 and σ(0) = ½ hold for the whole class. Own value at τ = 1.2: (h+h″)(10⁻⁶) = −0.382883 (c = 1) and −0.382889 (κ² = 9/7), against −0.382899.

---

## Own numerics (dps 60 unless noted)

**Perimeter, 8∫₀^{π/4}h du vs closed form − tail:**

| curve | τ | perimeter (50 digits) | difference | tail E |
|---|---|---|---|---|
| c = 1 | 0.3 | 8.371206555925532810555312805924362497395087116368 | 3.7·10⁻⁵⁸ | 0.0144782 |
| c = 1 | 0.9 | 11.34988672603641589852515129960334738298938003378 | 2.0·10⁻⁵⁶ | 14.374 (near the τ = 1 pole) |
| κ² = 9/7 | 0.3 | 8.354257372662867849879798894775448144815797846848 | 0 | 0.0024041 |
| κ² = 9/7 | 0.9 | 11.18730278932085423976972864298375489387023253798 | 0 | −0.099097 |

- The imaginary part of the closed form is 0 for κ² = 9/7.
- A full-circle quadrature over [0, 2π] with no symmetry reduction (dps 40, τ = 0.3) agrees to 1.8·10⁻⁴⁰ for both curves.

**min of h + h″ on u = (π/4)i/4000, i = 1…4000:**

| | τ = 0.3 | τ = 0.9 | τ = 1 | τ = 0.1 |
|---|---|---|---|---|
| c = 1 | 0.106946 (at π/4) | 0.209904 (u → 0, limit 0.207341) | 0.0027 (u → 0, limit 0) | 4.69·10⁻⁵ (at π/4) |
| κ² = 9/7 | 0.103880 (at π/4) | 0.209067 (u → 0) | 0.0019 (u → 0, limit 0) | 1.04·10⁻⁵ (at π/4) |

**Jumps of h′ at π/4:** 0.0433823 (c = 1, τ = 0.3), 0.781650 (c = 1, τ = 0.9), 0.0324483 and 1.011874 (κ² = 9/7).

**Curve from h.** I took 8 octants with N points each, using one-sided h′ at the diagonals so that the short edges are included.
- In all four cases the polygon has total turning exactly 1·2π and every normalised cross product is > 0. So it is a simple closed convex polygon.
- The length error falls as N⁻² (e.g. −2.1·10⁻⁵, −5.3·10⁻⁶, −1.3·10⁻⁶, −3.3·10⁻⁷ for N = 100…800).
- Double Richardson gives the perimeter to 10⁻¹⁶ (τ = 0.3) and 10⁻¹⁵ (τ = 0.9). For example, 8.37120655592553271 against 8.371206555925532811.

---

## What a paper may print

**As "proved":**
- Lemma 1, with the repaired step 2.
- Theorem A for 0 < τ < 1/Re a, and its continuation to all τ > 0 (pole cancellation). The c = ±1 case is stated as a limit, or by Theorem B.
- log R = Σ(a^k+b^k−1)ζ(k)τ^k/k, polynomiality in ab, and the table through τ⁷ (exact symbolic algebra).
- Theorem B with Λ = ½[q − W(e^q)], and the τ = 1 value.
- Admissibility ⟺ κ < κ_max. κ_max is the unique root of κe^{π/(2κ)} = 2 + κ², certified in (1.50781081090595335, 1.50781081090595336).
- The τ⁴ zero at κ² = 9/7, and the exact τ⁵ and τ⁶ coefficients.
- Theorem D (c₂…c₇ as printed), the ψ^{(2k)} coefficient 4^{1−k}η(2k+2), and the "iff" rigidity.
- The formula for h + h″; C² at the edge normals; the jump 8fσ/φ(V) > 0; and "∫h = perimeter", using the self-contained support-function argument above.
- E1 (log ψ convex, all real c, τ ∈ (0,1]).
- E2(i) for κ² = 9/7, with a computer-assisted (ball-arithmetic) certificate of one inequality. Say so.
- E2(ii) (τ ≤ 1/(1+κ tan κV); every member convex for τ ≤ 1/(1+δ)).
- E3 (no member is convex for τ > 1).

**Only as "numerical":**
- κ_s ≈ 1.18271 and κ_D ≈ 1.35542.
- The τ*(κ) table and the "min Dφ/M²" column (reword it: the infimum is 0 at u = 0).
- The grid minima of h + h″.
- The 15 quadrature agreements, and the decimal values of κ_max, V and φ(V), beyond the certified brackets above.

All scripts and logs are in this folder. Nothing outside it was modified.

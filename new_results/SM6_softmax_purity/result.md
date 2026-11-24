# SM6: are the tau^6 and tau^7 coefficients of the soft-max perimeter pure weight?

**Answer: no.** Both c6 and c7 contain one genuinely alternating Euler-sum constant that does not cancel, so neither is a rational combination of ordinary zeta products.

- c6 = (−173 ζ(6) + 44 ζ(3)² + 155 ζ(5) log 2 − 80 ζ(−5,−1)) / 32 = 0.042950827880284033929433804835135265510296727692058… (54 digits, Richardson)
- c7 = −(5/24) X + (8339/512) ζ(7) − (75/16) ζ(2)ζ(5) − (279/128) ζ(3)ζ(4) = 5.58574121431910295823282665831546806173150467798… (47 digits, Richardson)
  with X = ∫₀¹ (−log x)³ log(1−x) log(1+x) log(1−x²) dx/x, a weight-7 alternating constant.
- Exact reduction (PROVED, given Theorem 1 of LAY1 and the exact reduction below): c6 = 11/16 ζ(6) − ζ(3)²/2 + 40 γ₃₂ and c7 = 40 γ₃₃ + (3/4) ζ(3)ζ(4) + (2339/512) ζ(7). Lemma 5 and the η-moments do not remove γ₃₂ or γ₃₃. The relations are exactly symmetric under λ → 2μ−λ, and γ₃₂ and γ₃₃ are antisymmetric data.
- γ₃₂ = (1/16) Σ_{n≥2} (−1)^n H_{n−1}/n⁵ is a weight-6 linear alternating Euler sum. PSLQ puts it in the span of ζ(6), ζ(3)², ζ(5)log 2 and ζ(−5,−1), with ζ(−5,−1) coefficient ≠ 0.
- The p-norm ball (control) reduces completely, with no free γ, through c7 (checked at c2..c7).

Convention note. The task text writes h = (Σ exp(v·n/τ))^τ. The coefficients quoted (8, π²/3, 2ζ(3), …) and the pickled SM5 data belong to the LAY1 definition h = m(θ) [Σ_v exp((v·n − m)/τ)]^τ with m = |cos θ| + |sin θ|. I used that one throughout. The two differ by the factor e^{m} versus m, so only LAY1's has I = 8 + π²τ²/3 + ….

All files are in this folder. Scripts need mpmath 1.3 and sympy 1.14 (I used a throwaway venv here and deleted it).

---

## 1. Numerics (item 1)

`sm6_richardson.py` (c6), `sm6_richardson_c7.py` (c7). Input: the SM5 values I(τ) at τ = 1/(20·2^i), i = 0..17 (dps 170, quadrature error ≤ 9e-176). Subtract 8 + (π²/3)τ² + 2ζ(3)τ³ + (ζ(4)/4)τ⁴ + c5τ⁵ with c5 = (107/32)ζ(5) − ζ(2)ζ(3) (proved in SM5), divide by τ⁶:

Y(τ) = c6 + c7 τ + c8 τ² + …

Neville elimination with ratio 2. For c7 I apply Neville again to (Y − c6)/τ with c6 the Richardson value itself (so no use of any closed form).

Contamination check: the data contain exponentially small terms e^{−√2/τ·…}. They are < 1e−85 relative to the quantities used from n = 160 on. Starting at n = 20, n = 160 or n = 640 gives identical digits at every order shown below (the last k+1 points are what the entry uses), so the early points do no harm.

c6 (rows k = Richardson order, last entry of each column):

| k | c6 (first 50 digits) | change from order k−1 |
|---|---|---|
| 6  | 0.04295082788028403392943380483513531280190873098 | 1.3e-34 |
| 8  | 0.04295082788028403392943380483513526551030479020 | 1.6e-41 |
| 10 | 0.04295082788028403392943380483513526551029672772 | 6.0e-47 |
| 12 | 0.04295082788028403392943380483513526551029672769 | 4.7e-51 |
| 14 | 0.04295082788028403392943380483513526551029672769 | 7.4e-54 |
| 16 | 0.04295082788028403392943380483513526551029672769 | 2.6e-55 |
| 17 | 0.0429508278802840339294338048351352655102967276920578653499916 | 2.9e-55 |

The closed form of item 2 agrees with the k = 17 value to 1.5e-55. **c6 = 0.04295082788028403392943380483513526551029672769206 (claim 54 digits, noise floor ≈ 3e-55).**

c7:

| k | c7 | change from order k−1 |
|---|---|---|
| 6  | 5.5857412143191029582328266583154680618191500632 | 1.9e-30 |
| 8  | 5.5857412143191029582328266583154680617315047548 | 8.3e-38 |
| 10 | 5.5857412143191029582328266583154680617315046780 | 7.7e-44 |
| 12 | 5.5857412143191029582328266583154680617315046780 | 1.5e-48 |
| 14 | 5.5857412143191029582328266583154680617315046780 | 6.2e-52 |
| 16 | 5.5857412143191029582328266583154680617315046780 | 1.2e-53 |

The last stable digits are …046779809 (k = 13..17). The exact value of item 2 is …0467798015; the difference is 7.7e-49. The error of the Richardson c6 (1e-55) is amplified by 1/τ ≤ 2.6e6 and by Neville weights, so I claim **c7 = 5.585741214319102958232826658315468061731504678 (47 digits)**. A second run with c6 replaced by the 52-digit γ-formula value gives …046772751; it is worse by the truncation of c6 to 52 digits, so the honest figure is 47.

## 2. Exact reduction (item 2)

**Reduction, checked.** For 0 ≤ θ ≤ π/4 the terms e^{−2cos θ/τ}, e^{−2(cos θ+sin θ)/τ} are ≤ e^{−√2/τ}, so I = 8∫₀^{π/4}(cos θ + sin θ)(1+e^{−2 sin θ/τ})^τ dθ + O(e^{−√2/τ}). With u = sin θ: (cos θ+sin θ)dθ = (1 + u/√(1−u²))du, and ∫₀^{1/√2}(1+u/√(1−u²))du = 1 exactly, so
I/8 − 1 = ∫₀^{1/√2} (1+u/√(1−u²)) [(1+e^{−2u/τ})^τ − 1] du + O(e^{−√2/τ}).
Put u = τs and expand u/√(1−u²) = Σ_j C(2j,j)4^{−j}(τs)^{2j+1}. With G(λ,μ) = Σ_{a≥0,b≥1} γ_ab λ^a μ^b one gets (asymptotically, as a Watson-type expansion; the coefficient formula is what Theorem 1 of LAY1 gives, and the numerics below confirm it)

c_n/8 = γ_{0,n−1} + Σ_{j≥0} C(2j,j)4^{−j}(2j+1)! γ_{2j+1,n−2j−2}, with γ_{a,0} := 0,

that is,

- c4 = 8(γ03 + γ12)
- c5 = 8(γ04 + γ13 + 3γ31)
- **c6 = 8(γ05 + γ14 + 3γ32)**
- **c7 = 8(γ06 + γ15 + 3γ33 + 45γ51)**

(`sm6_verify.py`: the first two reproduce ζ(4)/4 and (107/32)ζ(5) − ζ(2)ζ(3) to 40 digits; the formulas for c6 and c7 reproduce the Richardson values to 1.4e-52 and 8e-49.)

**Relations.** Write ζ(2)=π²/6 etc. Lemma 5, in the form (after the Gamma expansion, with Q = Γ(1−λ/2)Γ(1+λ/2−μ)/Γ(1−μ), log Q = Σ_{k≥2} (ζ(k)/k)[(λ/2)^k + (μ−λ/2)^k − μ^k])

G(λ,μ) + G(2μ−λ,μ) = 2μ (Q − 1) / (λ(λ−2μ)),

is compared coefficient by coefficient in λ^aμ^b for all a+b ≤ 6 (this is all that differentiating in λ, μ can give), together with γ_{a1} = η(a+2)/2^{a+1} (η(n) = (1−2^{1−n})ζ(n)). `sm6_reduce.py` solves the linear system (sympy, π² and odd zetas as independent symbols). Solved in terms of the free ones (printed in `sm6_reduce_output.txt`; every one checked numerically at 80 digits, residuals ≤ 7e-80, in `sm6_verify_output.txt`):

| degree | solved | remains free |
|---|---|---|
| 3 | γ21 | γ12 (γ03 = −γ12 + π⁴/2880) |
| 4 | γ04 = −γ13 − π²ζ(3)/48 + 31ζ(5)/128, γ22, γ31 | **γ13** |
| 5 | γ05 = −γ14 + 2γ32 + 11π⁶/120960 − ζ(3)²/16, γ23 = −3γ32 − π⁶/24192 + ζ(3)²/32, γ41 | **γ14, γ32** |
| 6 | γ06 = −γ15 + 2γ33 + π⁴ζ(3)/960 − 31ζ(7)/256, γ24 = −3γ33 − …, γ42, γ51 | **γ15, γ33** |

Count: the relations are exactly the reflection-symmetric half, so the free data at degree d are the coefficients of odd powers of (λ−μ) in G_d (one fewer when the elementary γ_{d−1,1} falls on one).

**Cancellation.**
- c4 and c5: the free γ12 and γ13 cancel. This reproduces the proved ζ(4)/4 and the SM5 value of c5.
- c6 = 11/16 ζ(6) − ζ(3)²/2 + 40 γ32. γ14 cancels, γ32 does not (coefficient 40 = 8·(3+2)).
- c7 = 40 γ33 + (3/4)ζ(3)ζ(4) + (2339/512)ζ(7). γ15 cancels, γ33 does not.

So the first non-cancelling moment appears at the first order where an s³ moment meets g² (c6 = the a = 3, b ≥ 2 moments). Status: the algebra is PROVED (given Theorem 1 of LAY1 for the coefficient formula at N = 5, 6, whose (H2) I did not re-derive beyond N = 5 in SM5; the exponentially small tails are as in SM5 Section 2).

## 3. PSLQ (item 3)

`sm6_pslq6.py` (weight 6, dps 260), `sm6_pslq7.py`, `sm6_pslq7b.py` (weight 7, dps 420), `sm6_final.py`.

**Weight 6.**
- c6 (50 digits, Richardson) vs {ζ(6), ζ(3)²}: **no relation** (maxcoeff 1e6). Adding ζ(5)log 2: no relation.
- Basis of 13 (Fibonacci dimension of weight-6 alternating MZVs): ζ(6), ζ(3)², ζ(5)log2, ζ(4)log²2, ζ(3)log³2, ζ(2)log⁴2, log⁶2, Li6(1/2), Li5(1/2)log2, Li4(1/2)log²2, ζ(2)Li4(1/2), **ζ(2)ζ(3)log2** (missing from the list in the task; it is needed to make 13), and ζ(−5,−1) := Σ_{n>m≥1}(−1)^{n+m}/(n⁵m) = (1/24)∫₀¹(−log x)⁴ log(1−x)/(1+x) dx. These 13 are mutually independent at 200 digits (PSLQ, maxcoeff 1e9: none).
- c6 against these 13, from the Richardson value alone (50 digits, maxcoeff 1e3): **32 c6 + 173 ζ(6) − 44 ζ(3)² − 155 ζ(5)log2 + 80 ζ(−5,−1) = 0**, all other coefficients 0.
- The same relation from the exact c6 = 11/16 ζ(6) − ζ(3)²/2 + 40γ32 (γ32 by quadrature, dps 250): γ32 = (−39 ζ(6) + 12 ζ(3)² + 31 ζ(5)log2 − 16 ζ(−5,−1))/256, and the c6 relation above again, residual 1.5e-55 against Richardson. Also γ32 = (1/16)Σ_{n≥2}(−1)^n H_{n−1}/n⁵ numerically (20 digits).

So c6 = (−173 ζ(6) + 44 ζ(3)² + 155 ζ(5)log 2 − 80 ζ(−5,−1))/32. The alternating content is exactly ζ(5)log 2 and the depth-2 generator ζ(−5,−1); no Li_n(1/2) and no log^k 2 products appear.

**Weight 7.**
- c7 (47 digits) vs {ζ(7), ζ(2)ζ(5), ζ(3)ζ(4)}: **no relation** (maxcoeff 1e8). γ33 itself vs the same three at 300 digits: none (maxcoeff 1e12).
- Wider basis of 21 (Fibonacci dimension at weight 7): all 19 products of {log2, ζ2, ζ3, Li4(1/2), ζ5, Li5(1/2), Li6(1/2), ζ(−5,−1), ζ7, Li7(1/2)} of weight 7 plus two depth-3 integrals X1 = ∫(−log x)³ log²(1−x) log(1+x) dx/x and X2 = ∫(−log x)³ log(1−x) log²(1+x) dx/x. I tried all pairs of seven candidate integrals until the 21 numbers were independent (380 digits, maxcoeff 1e12; the first choices were reducible, see `sm6_pslq7_output.txt`).
- γ33 against these 21 plus γ33 (maxcoeff 1e14): **1536 γ33 + 8(X1+X2) − 450 ζ(7) + 180 ζ(2)ζ(5) + 45 ζ(2)²ζ(3) = 0**, all other coefficients 0. Hence γ33 = (−8 X + 450 ζ(7) − 180 ζ(2)ζ(5) − 45 ζ(2)²ζ(3))/1536 with X = X1 + X2 = ∫₀¹(−log x)³ log(1−x)log(1+x)log(1−x²) dx/x.
- Therefore c7 = −(5/24) X + (8339/512) ζ(7) − (75/16) ζ(2)ζ(5) − (279/128) ζ(3)ζ(4), verified: closed form − γ-formula = 1.7e-420, closed form − Richardson c7 = 7.7e-49.
- Directly from the Richardson c7: relation −1536 c7 + 25017 ζ(7) − 7200 ζ(2)ζ(5) − 3348 ζ(3)ζ(4) − 320 X = 0 (found at 47 digits, maxcoeff 1e5).
- X is independent of {ζ(7), ζ(2)ζ(5), ζ(3)ζ(4)} (300 digits, maxcoeff 1e12). Since the 21 numbers are independent, the representation of c7 is unique, so X cannot be absorbed.

The constant that enters at weight 7 is a depth-3 alternating Euler sum, not a product of lower-weight constants and not Li7(1/2) or log 2 powers: no Li_n(1/2) or log 2 term appears in c7 at all.

Caveat: "independent" and "not pure" are PSLQ statements, i.e. conditional on the standard expectation that these alternating constants are linearly independent over Q. I did not prove ζ(−5,−1) or X irrational or independent.

## 4. Control: the p-norm ball (item 4)

`sm6_squircle.py`, `sm6_squircle_num.py`.

**Exact form.** With φ = π/4 − θ, v = artanh(tan φ): dφ = dv/cosh 2v, cos φ + sin φ = √(1+tanh 2v), and m = cos θ on [0, π/4]. Then, with no truncation error at all,

I_p = 4√2 [1 + ∫₀^∞ w(v) ((1+e^{−2v/τ})^τ − 1) dv], w(v) = √(1+tanh 2v)/cosh 2v = 2√2 e^{−2v}(1+e^{−4v})^{−3/2}.

The soft-max has the same shape, I = 8 [1 + ∫₀^{1/√2} w_s(v)(…) dv] with w_s(v) = 1 + v/√(1−v²). Taylor coefficients of w: 1, 1, −5/2, −17/6, 121/24, 721/120, … (w_s: 1, 1, 0, 1/2, 0, 3/8, 0, 5/16). Then c_n = pref·Σ_k w_k k! γ_{k,n−1−k}.

**Reduction** with the same relations: every free γ cancels at every order tested, n = 2..7:

| n | p-norm c_n | free γ left |
|---|---|---|
| 2 | √2 ζ(2) | none |
| 3 | √2 ζ(3) | none |
| 4 | −(33√2/16) ζ(4) | none |
| 5 | −√2 (3/4 ζ(5) + 7/4 ζ(2)ζ(3)) | none |
| 6 | √2 (4039/256 ζ(6) − 7/8 ζ(3)²) | none |
| 7 | √2 (129/16 ζ(7) + 141/16 ζ(2)ζ(5) + 1071/64 ζ(3)ζ(4)) | none |

The n = 4, 5 values are the proved LAY1 and SM5 results. n = 6, 7 are new. Numerical check from the exact integral at τ = 1/(10·2^i), i ≤ 14, dps 150: Richardson c6 agrees with the formula to 3.6e-33, c7 to 1e-32 (the floor is the 150-digit quadrature noise amplified by 1/τ⁶; `sm6_squircle_num_output.txt`). Both are in the pure ordinary-zeta algebra by construction of the closed form.

**What cancels the free part there (my reading, not proved).** The p-norm support function (|cos θ|^q + |sin θ|^q)^{1/q} is real-analytic across the kink θ = π/4, so the two half-windows are analytic continuations of each other, and the folded half-line integral is half of a full-line Beta-type integral; this is the mechanism behind Lemma 5 (G(λ,μ) + G(2μ−λ,μ) is the full-line integral, up to the pole terms). The soft-max of LAY1, h = m[Σ e^{(v·n−m)/τ}]^τ, is only C² at its kink: the two branches (cos u + sin u)(1+e^{−2 sin u/τ})^τ and (cos u − sin u)(1+e^{2 sin u/τ})^τ differ by a factor 1 + O(u³). An s³ mismatch matches where the first free moment (a = 3, b ≥ 2) enters, at τ⁶. I have not turned this into a proof: it is consistent with the data (first impurity at c6, none before) but no test beyond that was run.

## 5. Status

| item | status |
|---|---|
| Exact reduction I/8 − 1 = Σ γ-expansion; c6 = 8(γ05+γ14+3γ32), c7 = 8(γ06+γ15+3γ33+45γ51) | **PROVED** (from Theorem 1 of LAY1 / exact change of variables; higher-order (H2) at N = 6, 7 not re-derived) |
| Relations from Lemma 5 (all degrees ≤ 6) and elementary moments; free set {γ12; γ13; γ14, γ32; γ15, γ33} | **PROVED** (sympy linear algebra; each relation checked ≤ 7e-80) |
| c6 = 11/16 ζ(6) − ζ(3)²/2 + 40γ32; c7 = 40γ33 + (3/4)ζ(3)ζ(4) + (2339/512)ζ(7) | **PROVED** (same caveat) |
| γ14 and γ15 cancel from c6 and c7; γ32 and γ33 do not cancel by any relation of Lemma 5 type | **PROVED** (they are antisymmetric data; the relations only fix the symmetric part) |
| c6 = 0.04295082788028403392943380483513526551029672769206 | **NUMERICAL** (54 digits) |
| c7 = 5.585741214319102958232826658315468061731504678 | **NUMERICAL** (47 digits) |
| c6 = (−173ζ(6) + 44ζ(3)² + 155ζ(5)log2 − 80ζ(−5,−1))/32 | **NUMERICAL** (PSLQ at 250 digits for γ32, relation also found on c6 itself at 50 digits; agreement 1.5e-55) |
| c7 = −(5/24)X + (8339/512)ζ(7) − (75/16)ζ(2)ζ(5) − (279/128)ζ(3)ζ(4) | **NUMERICAL** (PSLQ at 380 digits for γ33; agreement 7.7e-49) |
| c6, c7 are not pure ordinary-zeta | **NUMERICAL** (no relation with the pure bases; unique representation in an independent 13/21-element alternating basis, conditional on independence of ζ(−5,−1), X) |
| ζ(−5,−1), X are irrational / independent of ordinary MZVs | **OPEN** (standard expectation, not proved) |
| p-norm c6, c7 pure; free γ cancel through n = 7 | **PROVED** (reduction) + **NUMERICAL** (Richardson to 33 digits) |
| Structural reason (analyticity of h_p at the kink vs C² of the soft-max) | **OPEN** (heuristic) |
| Whether all c_n, n ≥ 6, of the soft-max contain an alternating constant | **OPEN** |
| c8 and beyond | **OPEN** (not attempted) |

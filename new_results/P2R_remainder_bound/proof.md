# Remainder of "Asymptotics at every angle": O(log n / n^3) is correct and sharp

Source: `paper/ieee/main.tex` l. 1399 (Theorem "Asymptotics at every angle"); notes
`research/07_push/P2_general_angle/claude/result.md` (T1-T3). Numerics: `check.py`, `check_output.txt` (this folder).

**Verdict.** The statement is true as written. The remainder is exactly of order log n / n^3:
E_n := L - [1+cos θ + 2z(μ/n)^ε + c_2/n^2] = -2 z a_2 (log n)/n^3 + O(z/n^3), with a_2 > 0 for every θ in (0,π),
so log n cannot be removed. The bound is uniform on the whole open interval: |E_n| ≤ C·sin²(θ/2)·log n/n^3 for all n ≥ 9, all θ in (0,π), with an absolute constant C.

## 0. Notation and what is taken from the notes

z = sin²(θ/2) in (0,1), c = cos θ = 1-2z, w = (1-z)/z, ℓ = log(1-z), u = wℓ, ε = 1/(n-1), n ≥ 9 (so ε ≤ 1/8).
Facts used (all in the notes):

* (F1) Exact formula (Thm 2, valid for n>1, n ≠ 1+1/j):
  L = 1 + c + R(ε) n^{-ε} + Σ_{m≥1} (-1)^m ε/(m-ε) Λ_m n^{-m}, with R(ε) = (πε/sin πε) Λ_ε(θ), Λ_ε = 2z·₂F₁(-ε,1+ε;2;z), Λ_m = (P_{m-1}(c)-P_{m+1}(c))/(2m+1).
* (F2) |Λ_m| ≤ 1-c = 2z for integers m ≥ 0 (from |P_m| ≤ 1 on [-1,1]); Λ_1 = 2z(1-z).
* (F3) R(ε) = ε M_θ(ε), M_θ(ε) = ∫_0^∞ k_c(s) s^{ε-1} ds, k_c(s) = sin²θ/(s+c+√(1+2cs+s²)) > 0 (Thm 1), 0<ε<1.
* (F4) g(ε) := R(ε)/(2z) = (πε/sin πε)·F(ε,z), F := ₂F₁(-ε,1+ε;2;z); g(0)=1. Put log g = a_1 ε + a_2 ε² + a_3 ε³ + ... (the notes' a_j), μ = e^{a_1}, c_2 = 2z(a_2 - 1 + z).
  The values a_1 = -1-u and a_2 = ζ(2)+1/2-Li_2(z)-u²/2 are Thm 3(a) of the notes (Prop. 4 there). Only these two values enter the statement. a_1 is checked below in two lines; a_2 is taken from the notes and re-checked numerically (§4).

## 1. Lemma A (uniform regularity of g in z, including z → 0 and z → 1)

**Lemma A.** There are constants B, g_min > 0 and K_3 < ∞, independent of z in [0,1], such that on ε in [0,1/8]: g ≥ g_min, and |(log g)'''(ε)| ≤ K_3.

*Proof.* Write F = 1 + Σ_{k≥1} t_k z^k, t_k = (-ε)_k(1+ε)_k/((2)_k k!), and (-ε)_k = -ε(1-ε)_{k-1}, (2)_k = (k+1)!.
For complex |ε| ≤ 1/4 put a = |ε|. Then Π_{j<k}(j+a) ≤ (k-1)! e^{a(1+log k)} and Π_{j≤k}(j+a) ≤ k! e^{a(1+log k)}, so
|t_k| ≤ a e^{2a} k^{2a} (k-1)! k!/((k+1)! k!) ≤ a e^{2a} k^{2a-2} ≤ (e^{1/2}/4) k^{-3/2}.
Since Σ k^{-3/2} = ζ(3/2) < 2.62, the series for F converges absolutely and uniformly for (ε,z) in {|ε|≤1/4}×[0,1], so F is holomorphic in ε, jointly continuous (also at z=1), and |F| ≤ 2.1.
The factor πε/sin πε is holomorphic and bounded on |ε| ≤ 1/4. Hence |g| ≤ B there. Cauchy's estimate on circles of radius 1/8 centred at ε_0 in [0,1/8] gives |∂_ε^j g| ≤ j! 8^j B, j ≤ 3, uniformly in z.

Positivity of g on [0,1/8]×[0,1]:
* z=0: F=1, g = πε/sin πε > 0.
* 0<z<1: for 0<ε<1/8, g = εM_θ(ε)/(2z) > 0 by (F3), and g(0)=1.
* z=1: Gauss' theorem (c-a-b = 1 > 0) gives F = 1/(Γ(2+ε)Γ(1-ε)), so g = Γ(1+ε)Γ(1-ε)/(Γ(2+ε)Γ(1-ε)) = 1/(1+ε) > 0 (this is the hairpin value of the notes).

g is continuous and strictly positive on the compact set [0,1/8]×[0,1], so g ≥ g_min > 0. Finally
(log g)''' = g'''/g - 3g'g''/g² + 2(g')³/g³ is bounded by a constant K_3 built from B and g_min. ∎

*Check of a_1.* t_k = -ε/(k(k+1)) + O(ε²), and Σ_{k≥1} z^k/(k(k+1)) = 1 + wℓ = 1+u. So F = 1 - ε(1+u) + O(ε²), and πε/sin πε = 1 + O(ε²), giving a_1 = -(1+u). Since ℓ < 0 and 0 < -ℓ(1-z)/z < 1 (because -log(1-z) < z/(1-z)), we have u in (-1,0), so **a_1 in (-1,0) and μ in (e^{-1},1)**.

**Lemma A'. a_2 in (0, ζ(2)+1/2).** a_2 = ζ(2) + 1/2 - Li_2(z) - u²/2 with 0 < Li_2(z) < ζ(2) and 0 < u² < 1 for z in (0,1). So 0 < a_2 < ζ(2)+1/2 =: A.

## 2. Lemma B (Taylor with explicit remainder)

For ε in [0,1/8], Taylor's theorem with Lagrange remainder applied to log g (which is C³ by Lemma A, with log g(0)=0) gives
log g(ε) = a_1 ε + a_2 ε² + r_3(ε), |r_3(ε)| ≤ K_3 ε³/6.
Put x = a_2 ε² + r_3. Then |x| ≤ (A + K_3/48) ε² =: A' ε², and e^x = 1 + x + ρ_2 with |ρ_2| ≤ x² e^{|x|}/2 ≤ (A'^2 e^{A'/64}/2) ε⁴. Therefore

  g(ε) = e^{a_1 ε} (1 + a_2 ε² + ρ(ε)),  ρ = r_3 + ρ_2,  |ρ(ε)| ≤ K_4 ε³  (K_4 depends only on K_3, A), uniformly in z in (0,1).

Consequently R n^{-ε} = 2z q (1 + a_2 ε² + ρ), with q := (μ/n)^ε = exp(-ε λ_n), λ_n := log n - a_1 in (log n, log n + 1).

## 3. The remainder bound (all n ≥ 9, all z in (0,1))

By (F1), (F4), Lemma B and Λ_1 = 2z(1-z), and since c_2/n² = 2z a_2/n² - 2z(1-z)/n²:

E_n = T_1 + T_2 + T_3 + T_4,
* T_1 = 2z a_2 (q ε² - n^{-2}) = 2z a_2 [(ε² - n^{-2}) - (1-q) ε²],
* T_2 = 2z q ρ(ε),
* T_3 = -Λ_1/(n(n-2)) + 2z(1-z)/n² (the m=1 term of the series, since ε/(1-ε) = 1/(n-2), minus its n^{-2} part),
* T_4 = Σ_{m≥2} (-1)^m ε/(m-ε) Λ_m n^{-m}.

Elementary facts for n ≥ 9: ε = (1/n)·n/(n-1) ≤ (9/8)/n; log n ≥ 2; ε² - n^{-2} = (2n-1)/(n²(n-1)²) ∈ [0, 3/n³]; 0 ≤ 1 - q ≤ ελ_n ≤ ε(log n + 1) ≤ 2ε log n (using 1-e^{-y} ≤ y and λ_n ≤ log n + 1 because a_1 > -1).

* |T_1| ≤ 2zA [3/n³ + 2ε³ log n] ≤ 2zA [3 + 2.85 log n]/n³ ≤ 9 A z (log n)/n³ (using log n ≥ 2).
* |T_2| ≤ 2z K_4 ε³ ≤ 3 K_4 z/n³ (q ≤ 1 since μ<1).
* |T_3| = 2z(1-z)·2/(n²(n-2)) ≤ 6z/n³.
* |T_4| ≤ 2z ε Σ_{m≥2} n^{-m} (by (F2) and m-ε ≥ 1) ≤ 2z ε n^{-2}/(1-1/n) ≤ 3z/n³.

Summing: **|E_n| ≤ z (9A + 3K_4 + 9) (log n)/n³ = C sin²(θ/2) (log n)/n³** for n ≥ 9 and θ in (0,π). ∎

So O(log n/n³) holds, with the constant uniform in θ.

### Sharpness (the true order)
Only T_1 carries a logarithm. (1-q)ε² = ε³λ_n + O(ε⁴λ_n²) = (log n - a_1)/n³ + O(log² n/n⁴), and T_2, T_3, T_4 and 2za_2(ε²-n^{-2}) are O(z/n³) without logs. Hence, uniformly in θ,

  **E_n = -2 z a_2 (log n)/n³ + O(z/n³)**, with -2za_2 < 0 by Lemma A'.

So E_n n³/log n → -2 sin²(θ/2)·a_2(θ) ≠ 0, and O(log n/n³) cannot be improved to O(1/n³). Carrying the n^{-3} terms through gives the full n^{-3} coefficient

  E_n n³ = d - 2 z a_2 log n + O(log² n/n),  d = 2z[a_2(2+a_1) + a_3] - z(3 - z - 2z²),

which is the notes' Q_3 (here: 2za_2·2 from ε²-n^{-2}; 2z(a_3 - a_2(log n - a_1)) from ε³ terms; -4z(1-z) from the m=1 term; Λ_2/2 = z(1-z)(1-2z) from m=2). This agrees with the table below.

## 4. Numerics (`check.py`, mpmath 1.4.1, 90 digits; full log in `check_output.txt`)

* The exact formula (F1) reproduces the literal arc length ∫_0^1 √(1+2c·n x^{n-1}+n²x^{2n-2}) dx to 1e-86 at θ = π/3, π/2, 2π/3, n = 3, 10.
* a_1, a_2, a_3 (closed forms of the notes) vs Cauchy-integral Taylor coefficients of log g on |ε| = 1/4: differences ≤ 1.3e-79.
* D(n) := E_n n³/log n; the prediction is -2za_2 + d/log n.

| θ | -2 z a_2 (limit) | n=1e2 | 1e3 | 1e4 | 1e5 | 1e6 |
|---|---|---|---|---|---|---|
| π/3 | -0.752429 | -0.59237 | -0.64431 | -0.67168 | -0.68792 | -0.69869 |
| π/2 | -1.322467 | -1.06069 | -1.14648 | -1.19121 | -1.21765 | -1.23515 |
| 2π/3 | -1.589546 | -1.29841 | -1.39561 | -1.44522 | -1.47434 | -1.49358 |

D(n) is bounded, monotone, and tracks the predicted -2za_2 + d/log n to ~1e-5 at n = 1e6 (prediction column in `check_output.txt`); it converges to the limit only at the 1/log n rate, as the d/log n term says. The next-order check r4 = (E n³ - (d - 2za_2 log n))·n/log² n stays bounded (≤ 0.54 for n ≤ 1e6), consistent with the O(log² n/n⁴) term. Extra angles θ = 0.1 and π-0.1 in the log show D/z of the same size, confirming the sin²(θ/2) scaling and uniformity near both ends.

## 5. Uniformity in θ, and the ends

* Compact subsets of (0,π): immediate. In fact the bound holds on all of (0,π) with C absolute (Lemma A is uniform on z in [0,1]).
* θ → 0 (z → 0): E_n = O(z log n/n³), i.e. the *relative* error is uniform. At θ = 0, L = 2 exactly and all explicit terms are 0 (z=0), consistent.
* θ → π (z → 1): the bound stays uniform. a_1 → -1, a_2 → 1/2, a_3 → -1/3 (from g = 1/(1+ε)), while u = wℓ → 0 only like (1-z)log(1-z), so the a_j are continuous but not smooth in θ at π. This does not matter, because Lemma A uses only z-uniform bounds on ε-derivatives. By continuity of L_θ(n) in θ, the inequality also holds at θ = π, where L_π(n) = 2(1-1/n)n^{-ε}.
* Only ε ≤ 1/8 (n ≥ 9) was used; for fixed bounded n the statement is vacuous as an O(·) claim.

## 6. Status table

| Step | Status |
|---|---|
| Exact formula (F1), Mellin transform (F3), \|Λ_m\| ≤ 2z | PROVED (notes T1, T2); formula re-verified against quadrature (1e-86) |
| Values a_1, a_2 (hence μ_θ, c_2) | a_1 PROVED here (§1 check); a_2 PROVED in notes (Prop. 4); NUMERICAL here (1e-79) |
| Lemma A: g analytic, bounded, positive, derivatives bounded uniformly in z in [0,1] | PROVED (this file; explicit M-test, Cauchy, Gauss at z=1) |
| Lemma A': 0 < a_2 < ζ(2)+1/2 | PROVED |
| Lemma B: Taylor with Lagrange remainder, exp expansion | PROVED |
| Bound on T_1..T_4, hence \|E_n\| ≤ C sin²(θ/2) log n/n³, n ≥ 9 | PROVED (constants explicit up to K_3, K_4, which are finite by Lemma A) |
| Sharpness E_n n³/log n → -2za_2 | PROVED (§3); NUMERICAL table consistent |
| n^{-3} coefficient d (= Q_3) | PROVED given a_3; NUMERICAL agreement ~1e-5 at n=1e6 |
| Uniformity in θ; behaviour at 0 and π | PROVED |

**Verdict: O(log n/n³) is PROVED as stated** (for each fixed θ in (0,π), uniformly in θ, with error ≤ C sin²(θ/2) log n/n³), and it is the exact order.

## 7. Suggested edits to the paper

1. The statement and the sketch are correct; no change to the claim. The sketch's use of "Taylor's theorem ... uniformly" needs Lemma A: the uniformity in θ comes from the z-uniform bound on ε-derivatives of g, which the sketch does not mention.
2. Optional strengthening of the statement: "uniformly in θ in (0,π), the error is at most C sin²(θ/2) log n/n³; it equals -2 sin²(θ/2) a_2 (log n)/n³ + O(1/n³) with a_2 > 0, so the logarithm is genuine." (a_2 > 0 follows from Li_2(z) < ζ(2), |u| < 1.)
3. In the sketch, "analytic and positive for 0 ≤ ε ≤ 1/2" is true but positivity should be cited from the Mellin integral of the positive function k_c (Thm 1), not from the ₂F₁ series (whose terms are negative for k ≥ 1).

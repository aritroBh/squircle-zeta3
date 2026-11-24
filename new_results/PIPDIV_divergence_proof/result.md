# PIPDIV — the 1/p-series of π_p diverges: proof and sharp large-order law

Notation throughout: h = 1/p, κ = h/(1−h), y₀ = iπ, Log is the principal logarithm (so Log y₀ = log π + iπ/2),
ψ = Γ'/Γ, γ = Euler's constant, (a)_m = a(a+1)…(a+m−1).
π_p/4 = 2^{−h} + Ê(h) with Ê = E/4 = hκ ∫₀^∞ e^{−hκs}(1+e^{−κs})^{−1−h}[(1+e^{−s})^h − 1] ds (the exact split of D12, Lemmas 1–2,
checked by the referee). π_p/4 ~ 1 + Σ_{k≥0} N_k h^{k+1}. Normalised coefficient: ν_k := N_k π^{k−1} / (2(k−2)!).

All constants hidden in O(·) below exist but are not computed (they come from compactness/maximum-modulus bounds).

---

## 1. Statements

### Theorem A (Laplace form and splitting of N_k) — PROVED

(i) For 0 < h < 1,
  Ê(h) = h Σ_{n≥1} C(h,n) ∫₀^∞ e^{−n y/h} H_n(y,h) dy,  H_n(y,h) = e^{(n+1)y} (1+e^y)^{−1−h},
where C(h,n) = binomial(h,n).

(ii) Ê has an asymptotic expansion to all orders, and
  N_k = (−log 2)^{k+1}/(k+1)! + T_k + U_k,
  T_k = Σ_{m=0}^{k−2} m! [y^m h^{k−2−m}] A(y,h),  A(y,h) := H_1(y,h) = e^{2y}(1+e^y)^{−1−h},
  |U_k| ≤ C (2^k + k (k−2)! 6^{−k}).
(T_k is the n = 1 term; U_k collects n ≥ 2.)

### Theorem B (complete large-order expansion) — PROVED

Let Φ(u,h) = e^{2u}((e^u−1)/u)^{−1−h} = Σ_l φ_l(h) u^l (so φ₀ = 1, φ₁ = (3−h)/2). For every L ≥ 1, as k → ∞,
  N_k = T_k^{(L)} + O_L( (k−2)! π^{−(k−1)} · k^{−L−1} log k ),
  T_k^{(L)} = 2 Re Σ_{l=0}^{L−1} Σ_{m=0}^{k−2} [h^{k−2−m}] ( (−1)^l φ_l(h) y₀^{l−1−m} e^{−h Log y₀} (1+h−l)_m ).
(Each T_k^{(L)} is an explicit finite sum of rational combinations of log π, π, γ and values ψ^{(j)}(integer).)

### Theorem C (two-term law; sharp form of the divergence) — PROVED

  N_k = 2 Re[ (k−2)! y₀^{1−k} + (k−3)! y₀^{2−k} (ψ(k−2) + γ − Log y₀) ] + O( (k−2)! π^{−(k−1)} log²k / k² ).

Equivalently, with σ_k = (−1)^{(k−1)/2} (k odd) and τ_k = (−1)^{(k−2)/2} (k even):

- **odd k:** N_k = σ_k · (2(k−2)!/π^{k−1}) · [ 1 + π²/(2(k−2)) + O(log²k/k²) ];
- **even k:** N_k = τ_k · (2(k−3)!/π^{k−2}) · [ ψ(k−2) + γ − log π + O(log²k/k) ],
  i.e. ν_k = τ_k · π(log k + γ − log π)/k · (1 + O(log k/k)).

The sign of N_k (for large k) has period four: + + − − for k ≡ 1, 2, 3, 0 (mod 4).

**Three-term refinement — PROVED** (same argument with L = 2 and one more term of each sum; X_j := ψ(j) + γ − log π):
- odd k: ν_k = σ_k [ 1 + π²/(2(k−2)) − π²( X_{k−3}² − 5π²/12 − 3 + ψ'(k−3) ) / (2(k−2)(k−3)) + O(log³k/k³) ];
- even k: ν_k = τ_k · (π/(k−2)) [ X_{k−2} + π² X_{k−3}/(2(k−3)) + O(log³k/k²) ].

### Corollary D (radius zero) — PROVED

|N_k|^{1/k} = (k/(eπ))(1+o(1)) as k → ∞ (both parities). Hence Σ N_k h^{k+1} has radius of convergence 0.
The same holds at the other end p → 1: with η = 1/q − 1, π_q/4 = 1 + Σ_k (−1)^k N_{k−1} η^k (duality, D12 §1.5), so that series also has radius 0.

### Remark (Borel picture) — heuristic, not used

The n-th term of Theorem A(i) is a Laplace transform in the variable t = ny with kernel e^{−t/h}; its Borel singularities sit at
t = ±iπ n(2ν+1). The nearest ones, ±iπ (n = 1), give the (k−2)!/π^{k−1} growth; their residue is e^{2y₀}/e^{y₀} = e^{±iπ} = −1,
which fixes the sign; and the exponent −1−h turns each pole into a branch point (y₀−y)^{−1−h}, which produces the
log k in the even-k law. Optimal truncation error is then of size e^{−π/h} = e^{−πp} (consistent with D12's p = 10 test).

---

## 2. Proofs

### 2.1 Lemma 1 (Laplace form; Theorem A(i))
For 0 < h < 1 and s > 0, (1+e^{−s})^h − 1 = Σ_{n≥1} C(h,n) e^{−ns}. For 0<h<1 the signs of C(h,n) are (−1)^{n−1}, so
Σ_{n≥1}|C(h,n)| = −Σ_{n≥1} C(h,n)(−1)^n = 1 − (1−1)^h = 1. The integrand of Ê is therefore dominated by
hκ e^{−hκs}·1, which is integrable, and Fubini gives Ê = hκ Σ_n C(h,n) ∫₀^∞ e^{−hκs}(1+e^{−κs})^{−1−h} e^{−ns} ds.
Substitute y = κs. Since 1/κ = 1/h − 1, e^{−ns} = e^{−ny/h} e^{ny}, and
e^{−hy}(1+e^{−y})^{−1−h} e^{ny} = e^{(n−h)y} e^{(1+h)y}(1+e^y)^{−1−h} = H_n(y,h). □
(Checked numerically to 35 digits at p = 4, 10 and against π₁₀: pd3 (a).)

### 2.2 Lemma 2 (asymptotic expansion; Theorem A(ii))
*Analytic input.* On |w| < π, 1+e^w ≠ 0; let ℓ(w) be the branch of log(1+e^w) with ℓ(0) = log 2 (real on the real axis).
Put C₀ = max{|exp(−(1+h)ℓ(w))| : |w| ≤ 3, |h| ≤ 1/2} < ∞. For n ≥ 1, |t| ≤ 3 we have |t/n| ≤ 3, so
|H_n(t/n,h)| ≤ e^{6}C₀ =: M₀; for n ≥ 2, |t| ≤ 6 we have |t/n| ≤ 3 and |H_n(t/n,h)| ≤ e^9C₀ =: M₁.
Let b_{n,m}(h) = [t^m] H_n(t/n,h). Cauchy's estimates give, for |h| ≤ 1/2,
  |b_{n,m}(h)| ≤ M₀ 3^{−m} (n ≥ 1),  |b_{n,m}(h)| ≤ M₁ 6^{−m} (n ≥ 2).        (2.1)
Each b_{n,m} is analytic in |h| < 1/2 (Cauchy integral in t of a function analytic in h).
*Binomial bound.* |C(h,n)| = (|h|/n) ∏_{i<n}|1 − h/i| ≤ (|h|/n) e^{|h|H_{n−1}} ≤ |h| e^{|h|} n^{|h|−1}; for |h| ≤ 1/2,
|C(h,n)| ≤ e^{1/2}|h| n^{−1/2}, so Σ_n |C(h,n)|/n ≤ e^{1/2}ζ(3/2)|h|.                      (2.2)
*Watson step, uniform in n.* Let 0 < h ≤ 1/4 be real and M ≥ 1. Write n I_n(h) = ∫₀^∞ e^{−t/h} H_n(t/n,h) dt (t = ny).
On [0,2]: |H_n(t/n,h) − Σ_{m<M} b_{n,m} t^m| ≤ Σ_{m≥M} M₀(t/3)^m ≤ 3M₀(t/3)^M by (2.1).
On [2,∞): for real h ≥ 0, H_n(y,h) ≤ e^{(n+1)y}e^{−(1+h)y} ≤ e^{ny}, so ∫₂^∞ e^{−t/h}H_n(t/n,h)dt ≤ ∫₂^∞ e^{−t/(2h)}dt = 2h e^{−1/h};
also ∫₂^∞ e^{−t/h} t^m dt ≤ e^{−1/h} m!(2h)^{m+1}. Together:
  n I_n(h) = Σ_{m<M} m! h^{m+1} b_{n,m}(h) + ε_n(h),  |ε_n(h)| ≤ C_M h^{M+1}, uniformly in n.
Multiply by hC(h,n)/n and sum (Σ|C(h,n)| = 1 for real h):
  Ê(h) = Σ_{m<M} m! h^{m+2} D_m(h) + O(h^{M+2}),  D_m(h) := Σ_{n≥1} C(h,n) b_{n,m}(h)/n.
By (2.1)–(2.2) the series for D_m converges uniformly on |h| ≤ 1/2, so D_m is analytic there with |D_m(h)| ≤ C 3^{−m}|h|;
in particular D_m(0) = 0. Expanding each D_m in its Taylor series (radius ≥ 1/2) gives the asymptotic expansion of Ê to every
order. By uniqueness of asymptotic expansions, with d_{m,r} = [h^r]D_m,
  N_k = (−log2)^{k+1}/(k+1)! + Σ_{m=0}^{k−2} m! d_{m,k−1−m}.                       (2.3)
*Split n = 1 / n ≥ 2.* C(h,1) = h and b_{1,m}(h) = [y^m]H_1(y,h) =: a_m(h), so the n = 1 part of D_m is h a_m(h) and contributes
T_k = Σ_m m![h^{k−2−m}]a_m(h). The rest D_m^{≥2} satisfies |D_m^{≥2}(h)| ≤ C₁6^{−m} on |h| ≤ 1/2 by (2.1)–(2.2), hence
|[h^r]D_m^{≥2}| ≤ C₁6^{−m}2^r and
|U_k| ≤ C₁ Σ_{m=0}^{k−2} m! 6^{−m} 2^{k−1−m} = C₁2^{k−1} Σ_{m≤k−2} m!12^{−m} ≤ C₁ 2^{k−1}(13 + (k−2)(k−2)!12^{−(k−2)})
(m!12^{−m} ≤ 1 for m ≤ 12 and is increasing for m ≥ 11), i.e. |U_k| ≤ C₁(13·2^{k−1} + 72 k (k−2)! 6^{−k}). □
Relative to (k−2)!π^{−(k−1)} this is O(k(π/6)^k) → 0. (Numerically U_k ≈ ±(k−2)!/(2π)^{k−1}, from the n = 2 singularities at
t = ±2πi: pd2 Table 1.) The term (log2)^{k+1}/(k+1)! is negligible too.

### 2.3 Local structure of A at ±iπ
Let D = {|y| < 3π/2} minus the two slits σ_± = {±iπ·t : t ≥ 1}. D is simply connected, contains 0, and 1+e^y ≠ 0 on D
(zeros are at odd multiples of iπ). So ℓ(y) = log(1+e^y) (ℓ(0) = log 2) and A(y,h) = e^{2y}exp(−(1+h)ℓ(y)) are analytic on D,
for every complex h.

Let φ(u) = (e^u−1)/u (zeros at u = 2πiν, ν ≠ 0), log φ the branch on |u| < 2π with log φ(0) = 0, and
Φ(u,h) = e^{2u}exp(−(1+h)log φ(u)), analytic on |u| < 2π. For y₀ ∈ {iπ, −iπ} define on C∖y₀[1,∞)
  (y₀−y)^{a} := exp(a Log y₀) · exp(a Log(1−y/y₀)).

**Lemma 3.** On V_± := {|y ∓ iπ| < π/2}∖σ_±, A(y,h) = (y₀−y)^{−1−h} Φ(y−y₀, h) with y₀ = ±iπ.
*Proof.* Both sides are analytic on V_± (V_± ⊂ D; the disc minus a radius is connected). Near y = 0 both are exp of (2y − (1+h)·[a
logarithm of 1+e^y]), since 1+e^y = (y₀−y)φ(y−y₀) when e^{y₀} = −1. It suffices that the two logarithms agree at one point of a
path into V_±: on the segment y = iθ (y₀ = iπ, 0 ≤ θ < π), Log y₀ + Log(1−y/y₀) = log(π−θ) + iπ/2 and log φ(iθ−iπ) has argument
−(π−θ)/2 ∈ (−π/2, 0] (φ(−iθ') = e^{−iθ'/2}·2 sin(θ'/2)/θ'), with no branch crossing; at θ = 0 the sum is
log π + iπ/2 + log(2/π) − iπ/2 = log 2 = ℓ(0). The case y₀ = −iπ is the complex conjugate. Identity theorem. □
(Checked numerically at complex h and several y: pd3 (b), agreement 10^{−31}.)

Define S_l^±(y,h) := φ_l(h)(y−y₀)^l (y₀−y)^{−1−h} = (−1)^l φ_l(h) (y₀−y)^{l−1−h}, analytic on C∖y₀[1,∞), and
R_L := A − Σ_± Σ_{l<L} S_l^±. Since [x^m](1−x)^{−β} = (β)_m/m!,
  m! [y^m] S_l^± = (−1)^l φ_l(h) y₀^{l−1−m} e^{−h Log y₀} (1+h−l)_m.                    (2.4)

**Lemma 4 (Hankel bound).** For each L ≥ 1 there is C_L with, for all m ≥ 0 and |h| ≤ 1/4,
  |[y^m] R_L(·,h)| ≤ C_L max(m,1)^{|h|−L} π^{−m},  and  |[y^m] R_L(·,0)| ≤ C (2π)^{−m}.
*Proof.* (h = 0.) A(y,0) = e^{2y}/(1+e^y) is meromorphic with simple poles in |y| < 3π only at ±iπ, residue e^{2y₀}/e^{y₀} = −1;
S₀^±(y,0) = 1/(y₀−y) has residue −1; S_l^±(y,0) (l ≥ 1) is a polynomial. So R_L(·,0) is analytic on |y| < 3π; Cauchy on |y| = 2π.

(General h.) *Local bound.* M_Φ := max{|Φ(u,h)| : |u| = π, |h| ≤ 1/4} < ∞, so |φ_l(h)| ≤ M_Φπ^{−l} and
|Φ − Σ_{l<L}φ_l u^l| ≤ 2M_Φ(|u|/π)^L for |u| ≤ π/2. Also |(y₀−y)^{−1−h}| ≤ C|u|^{−1−Re h} (|e^{−h Log w}| ≤ |w|^{−Re h}e^{π|Im h|}).
Hence on V_+ (y₀ = iπ), with P_+ := A − Σ_{l<L}S_l^+,
  |P_+(y,h)| ≤ C|u|^{α},  α := L − 1 − Re h ∈ [L − 5/4, L − 3/4] ⊂ (−1, ∞),  u = y − iπ,
and R_L = P_+ + g_+ with g_+ := −Σ_{l<L}S_l^− analytic on the whole disc |y − iπ| < π/2 (its singularities lie on σ_−).
Same at −iπ.
*Contour.* Let m ≥ m₀(L) (m₀ ≥ 4, m₀ > 2L+2). [y^m]R_L = (1/2πi)∮_Γ R_L y^{−m−1}dy, where Γ is the circle |y| = 5π/4 except
that at ±5iπ/4 it runs in along the left side of the slit to y₀(1+1/m), around the circle |y−y₀| = π/m, and back out along the
right side (Cauchy's theorem on D, letting the contour approach the slits; R_L has continuous boundary values on each side of
σ_± away from y₀). Let Γ_+ be the part of Γ inside the disc |y−iπ| < π/2.
- Arcs of |y| = 5π/4 outside both discs: R_L is analytic near this compact set and bounded uniformly in |h| ≤ 1/4; contribution
  O((5π/4)^{−m}).
- ∫_{Γ_+} g_+ y^{−m−1}dy equals the integral over the arc of |y| = 5π/4 joining the endpoints of Γ_+ straight across the slit
  (g_+ y^{−m−1} is analytic in the disc, which avoids 0): O((5π/4)^{−m}).
- ∫_{Γ_+} P_+ y^{−m−1}dy: on the two outer arcs inside the disc |u| ≥ π/4, contribution O((5π/4)^{−m}); on the small circle
  |P_+| ≤ C(π/m)^α, |y|^{−m−1} ≤ (π(1−1/m))^{−m−1} ≤ e^{5/3}π^{−m−1}, length 2π²/m: O(m^{−α−1}π^{−m}); on each slit side
  y = iπ(1+t), 1/m ≤ t ≤ 1/4: ≤ Cπ^{α−m}∫₀^∞ t^α(1+t)^{−m−1}dt = Cπ^{α−m}Γ(α+1)Γ(m−α)/Γ(m+1) = O(m^{−α−1}π^{−m}),
  uniformly for α in the compact range above.
Since −α−1 = Re h − L ≤ |h| − L, this gives the bound for m ≥ m₀; for m < m₀ use Cauchy on |y| = 1. □

**Proposition 5.** T_k = T_k^{(L)} + O_L((k−2)! π^{−(k−1)} k^{−L−1} log k).
*Proof.* By (2.4), T_k − T_k^{(L)} = Σ_{m=0}^{k−2} m! [h^{i}] [y^m]R_L(·,h) with i = k−2−m ([y^m]R_L is analytic in |h| ≤ 1/4).
i = 0: by Lemma 4, ≤ C(k−2)!(2π)^{−(k−2)}. i ≥ 1: Cauchy in h on |h| = ρ := 1/log k (k ≥ e⁴), and m^{ρ} ≤ k^{ρ} = e:
  |m![h^i][y^m]R_L| ≤ e C_L m! (log k)^i max(m,1)^{−L} π^{−m} =: t_i.
For m ≥ m₁ := ⌈2^{L+1}π log k⌉ the ratio t_{i+1}/t_i ≤ (π log k/m)·2^L ≤ 1/2, so these terms sum to ≤ 2t₁ =
2eC_L(k−3)!(log k)(k−3)^{−L}π^{−(k−3)}. The ≤ m₁ terms with m < m₁ are each ≤ C m₁!(log k)^{k}, and
m₁·m₁!(log k)^k = o((k−3)!π^{−k}k^{−L−1}) since log (k−3)! ~ k log k ≫ k log log k. Collect, using (k−3)! = (k−2)!/(k−2). □

**Proof of Theorem B.** (2.3), Lemma 2 (U_k and the log 2 term are o of any power of k times (k−2)!π^{−(k−1)}), Proposition 5. □

### 2.4 Proof of Theorem C
*Gamma bound.* For M ≥ 1 and |h| ≤ 1/4: |Γ(M+h)| ≤ Γ(M+Re h) (|Γ(x+iy)| ≤ Γ(x) for x > 0, from Euler's integral). For
0 ≤ σ ≤ 1/4, log-convexity gives Γ(M+σ) ≤ Γ(M)^{1−σ}Γ(M+1)^σ = Γ(M)M^σ. For −1/4 ≤ σ < 0 and M ≥ 2, Γ(M+σ) ≤
Γ(M−1)^{|σ|}Γ(M)^{1−|σ|} = Γ(M)(M−1)^{−|σ|} ≤ Γ(M); for 1 ≤ M < 2, Γ(M+σ)/Γ(M) ≤ Γ(3/4)/min Γ < 1.39. So
|Γ(M+h)/Γ(M)| ≤ 1.4 M^{|h|}. (Grid check: max ratio 1.2254·M^{1/4}, pd3 (c).) Also |1/Γ(1+h)| ≤ C_Γ and
|e^{−h Log y₀}| ≤ e^{|h|·|Log y₀|} ≤ e^{0.5} on |h| ≤ 1/4.

Take L = 1: φ₀ = 1 and (1+h)_m = Γ(m+1+h)/Γ(1+h), so with i = k−2−m, M = k−1−i,
  T_k^{(1)} = 2 Re Σ_{i=0}^{k−2} y₀^{−(k−1−i)} [h^i] F_M(h),  F_M(h) := e^{−h Log y₀} Γ(M+h)/Γ(1+h).
[h^0]F_M = Γ(M); [h^1]F_M = Γ(M)(ψ(M) + γ − Log y₀) (since (log Γ(1+h))' = −γ at 0). The i = 0, 1 terms are exactly the
displayed two-term expression. For i ≥ 2, Cauchy on |h| = 1/log k with the Gamma bound gives
|y₀^{−(k−1−i)}[h^i]F_M| ≤ C(log k)^i (k−2−i)! π^{−(k−1−i)}; summing exactly as in Proposition 5 gives
O((k−4)! log²k π^{−(k−3)}) = O((k−2)!π^{−(k−1)} log²k/k²). Add Theorem B with L = 1 (error O(k^{−2}log k) relative). □

*Parities.* Put y₀ = iπ and X = ψ(k−2) + γ − log π, so ψ(k−2) + γ − Log y₀ = X − iπ/2. Then
2Re[(k−2)!y₀^{1−k}] = 2(k−2)!π^{1−k} Re(i^{1−k}) and 2Re[(k−3)!y₀^{2−k}(X − iπ/2)] = 2(k−3)!π^{2−k}[Re(i^{2−k})X + (π/2)Im(i^{2−k})].
Odd k: i^{1−k} = σ_k, i^{2−k} = iσ_k, so the bracket is σ_kπ/2 and N_k = σ_k(2(k−2)!/π^{k−1})(1 + π²/(2(k−2))) + error.
Even k: Re(i^{1−k}) = 0, i^{2−k} = τ_k, so N_k = τ_k (2(k−3)!/π^{k−2}) X + error, and the error is O(log²k/k) relative to
2(k−3)!/π^{k−2}. □

*Three-term refinement.* Use L = 2 in Theorem B (error O(k^{−3}log k) relative). Add the i = 2 term of the l = 0 sum,
[h²]F_M = Γ(M)·½[(ψ(M)+γ−Log y₀)² + ψ'(M) − ζ(2)] with M = k−3, and the i = 1 term of the l = 1 sum,
−φ₁(0)·y₀^{−(k−3)}(k−4)! = −(3/2)(k−4)! y₀^{−(k−3)} (from (h)_m = h(m−1)! + O(h²), m = k−3). The remaining terms
(i ≥ 3 for l = 0, i ≥ 2 for l = 1) are bounded by the same Cauchy argument: O(log³k/k³) and O(log k/k³) relative. Taking real
parts as above (Re(X − iπ/2)² = X² − π²/4, Im(X − iπ/2)² = −πX) gives the displayed odd and even formulas. □

### 2.5 Proof of Corollary D
By Theorem C, for odd k, |N_k| ≥ (k−2)!/π^{k−1} for large k; for even k, |N_k| ≍ (k−3)! log k/π^{k−2}. Stirling gives
|N_k|^{1/k} = (k/(eπ))(1+o(1)) in both cases, so limsup |N_k|^{1/k} = ∞ and the radius is 0. The p → 1 statement follows because
the coefficients there are ±N_{k−1}. □

---

## 3. Numerical checks

N_k computed independently by `pd1_Nk.py` (J-integral table + exact y-series of e^{−hy}(1+e^{−y})^{−1−h}; not the referee's
node-wise series method). Two runs with different precision, quadrature degree and cut-off (110 dps/96 nodes/S=760 and
90 dps/48 nodes/S=700) agree to 9.5×10^{−80} relative for all k ≤ 100; N_0..N_5 match the closed forms to 10^{−110}; N_21, N_30
match the referee's r5 to its 13 printed digits. Lemma 2 at order 6 (direct π_p at p = 200, 400, 800) converges to N_5: pd3 (d).

**T_k/N_k** (the n = 1 term carries everything): 0.99934 (k=13), 0.99999923 (k=21), 1 − 1.6×10^{−12} (k=40), 1 to 12+ digits beyond.
U_k(2π)^{k−1}/(k−2)! stays bounded (±2.06 at k = 21, ±1.17 at k = 61, −0.215 at k = 100), as the n = 2 singularities predict.

**Table.** ν_k = N_kπ^{k−1}/(2(k−2)!); "2-term", "3-term" are the laws of Theorem C; r₂ = (ν_k − ν₂)k²/log²k and
r₃ = (ν_k − ν₃)k³/log³k should stay bounded (they do); the last four columns are |T_k − T_k^{(L)}|/(2(k−2)!/π^{k−1}) for
L = 1..4 (Theorem B predicts O(k^{−L−1}log k)). Planner's fits for comparison: odd 1 + 5/k, even 11/k.

| k | ν_k | 2-term | 3-term | r₂ | r₃ | L=1 | L=2 | L=3 | L=4 | planner |
|---|---|---|---|---|---|---|---|---|---|---|
| 21 | 1.2833115 | 1.2597264 | 1.2855403 | 1.12 | −0.73 | 5.4e−2 | 2.3e−3 | 1.6e−3 | 4.7e−5 | 1.238 |
| 22 | 0.5086478 | 0.3774639 | 0.4733540 | 6.65 | 12.72 | 2.5e−2 | 6.8e−3 | 1.2e−4 | 2.6e−4 | 0.500 |
| 25 | 1.2187269 | 1.2145566 | 1.2224824 | 0.25 | −1.76 | 3.5e−2 | 1.1e−3 | 6.3e−4 | 5.9e−6 | 1.200 |
| 26 | 0.4298209 | 0.3389728 | 0.4104806 | 5.79 | 9.83 | 1.4e−2 | 3.6e−3 | 5.7e−5 | 8.3e−5 | 0.423 |
| 29 | 1.1787854 | 1.1827705 | 1.1823323 | −0.30 | −2.27 | 2.4e−2 | 5.7e−4 | 3.0e−4 | 6.1e−7 | 1.172 |
| 30 | 0.3753849 | 0.3081820 | 0.3637491 | 5.23 | 7.98 | 8.8e−3 | 2.1e−3 | 2.9e−5 | 3.4e−5 | 0.367 |
| 39 | −1.1239298 | −1.1333730 | −1.1262296 | 1.07 | 2.77 | 1.2e−2 | 1.7e−4 | 7.2e−5 | 1.9e−7 | 1.128 |
| 40 | −0.2903390 | −0.2527210 | −0.2861291 | −4.42 | −5.37 | 3.6e−3 | 7.6e−4 | 7.4e−6 | 6.0e−6 | 0.275 |
| 59 | −1.0781631 | −1.0865755 | −1.0791540 | 1.76 | 3.00 | 5.0e−3 | 3.2e−5 | 1.1e−5 | 3.1e−8 | 1.085 |
| 60 | −0.2059596 | −0.1887276 | −0.2049845 | −3.70 | −3.07 | 1.1e−3 | 2.0e−4 | 1.0e−6 | 6.1e−7 | 0.183 |
| 79 | −1.0576617 | −1.0640883 | −1.0581744 | 2.10 | 3.03 | 2.7e−3 | 1.0e−5 | 3.2e−6 | 6.6e−9 | 1.063 |
| 80 | −0.1624135 | −0.1523580 | −0.1620889 | −3.35 | −1.98 | 4.5e−4 | 7.8e−5 | 2.6e−7 | 1.3e−7 | 0.138 |
| 99 | −1.0458923 | −1.0508742 | −1.0461943 | 2.31 | 3.02 | 1.7e−3 | 4.2e−6 | 1.2e−6 | 1.9e−9 | 1.051 |
| 100 | −0.1352790 | −0.1286239 | −0.1351508 | −3.14 | −1.31 | 2.4e−4 | 3.9e−5 | 8.9e−8 | 3.9e−8 | 0.110 |

Reading: the L = 1 column is ≈ 16/k² at odd k (k = 59: 5.0e−3·59² = 17; k = 99: 16) and the higher L columns fall like higher
powers of 1/k, as Theorem B says. The planner's even-k fit 11/k drifts away from the data (0.110 vs 0.135 at k = 100, because the true law has a log k); its
odd-k fit 1 + 5/k is the first-order law 1 + π²/(2(k−2)) (π²/2 = 4.93). Full tables: `pd2_log_K100.txt`.

---

## 4. Statement for the paper (first person plural)

We prove that the 1/p-expansion of π_p diverges for every p: its coefficients satisfy
N_k = (−1)^{(k−1)/2}·2(k−2)!/π^{k−1}·(1 + π²/(2(k−2)) + O(log²k/k²)) for odd k, while for even k they are smaller by a factor
π(log k + γ − log π)/k, with signs of period four. The growth comes from the poles of 1/(1+e^{y}) at y = ±iπ, which become
branch points because of the exponent −1−1/p; in particular the series has radius of convergence 0.

---

## 5. Scripts (all in this folder; run with `../../../.venv/bin/python`)

- `pd1_Nk.py K DPS DEG S` — N_0..N_K from the exact formula N_k = (−log2)^{k+1}/(k+1)! + Σ_j j![h^k](κ^{j+1}g_j B_j)
  (J_{r,j} by composite Gauss–Legendre; g_j by exact series algebra). Writes `pd1_N_K*_*.txt`; logs `pd1_log_K60.txt`,
  `pd1_log_K100a.txt`, `pd1_log_K100b.txt` (two independent settings, K = 100).
- `pd2_asymptotics.py NFILE K` — computes T_k (n = 1 term) exactly from the series of e^{2y}(1+e^y)^{−1−h}, U_k = rest,
  the singular approximants T_k^{(L)} (L = 1..4, φ_l by series algebra, Pochhammer polynomials in h), and the two- and three-term
  laws; prints Tables 1–2. Logs `pd2_log_K60.txt`, `pd2_log_K100.txt`.
- `pd3_identities.py [NFILE]` — checks of the proof's identities: (a) Laplace form of Lemma 1 vs the s-integral and π₁₀;
  (b) the branch identity of Lemma 3 at complex h; (c) the Gamma-ratio bound; (d) Lemma 2's remainder at order 6. Log `pd3_log.txt`.

## 6. Where the planner's heuristic was wrong

1. **Sign.** Expanding only 1/(1+e^{−y}) (y = κs) against η(n+2) gives N_k ≈ c_{k−2}(k−2)!η(k), whose sign (−1)^{(k−3)/2} is the
   *opposite* of the true one. The corrections are not small: κ^{j+1} = h^{j+1}(1−h)^{−j−1}, and summing the (1−h)^{−j−1}
   factors over j multiplies the singular part by e^{y₀} = e^{±iπ} = −1. Equivalently, the right Borel variable is 1/h, not 1/κ:
   e^{−ny/κ} = e^{−ny/h}e^{ny}, and the n = 1 symbol is e^{2y}(1+e^y)^{−1−h}, with residue −1 at ±iπ. The magnitude
   2(k−2)!/π^{k−1} was right.
2. **Even k.** The suppression is not 11/k. It is π(log k + γ − log π)/k (times 1 + O(log k/k)): the exponent −1−h makes the
   poles branch points (y₀−y)^{−1−h}, and the first h-derivative brings log k. The even-k leading term comes entirely from that
   log, not from the other corrections the planner listed.
3. **Odd-k correction.** The O(1/k) term is exactly π²/(2(k−2)) (the planner's "5/k"); the error after it is O(log²k/k²), so the
   error term after the leading law is O(1/k), not O(log k/k).
4. **The corrections the planner listed** (exponent −1−h, e^{−hκs}, κ = h+h²+…, higher terms of (1+e^{−s})^h − 1) are all
   handled exactly by Theorem A: the first three are absorbed into one symbol e^{(n+1)y}(1+e^y)^{−1−h}; the last is the n ≥ 2
   sum, which only contributes O((π/6)^k) relative (numerically (1/2)^k).

---

## Planner note, 2026-11-22 (after the referee pass `../PIPDIVR_referee_divergence/REFEREE.md`)

Verdict: Theorems A, B, C and Corollary D (p → ∞) PROVED; two small repairs, both accepted.

- **R1 (Lemma 3).** The two logarithms are compared on W = D ∩ {|y − y₀| < 2π} ∖ σ_±, which is connected and contains both 0 and V_±. On W two continuous logarithms of 1 + e^y differ by a constant in 2πiℤ, and the constant is 0 at y = 0. (The original text said "near y = 0", but 0 is not in V_±.)
- **R2 (Corollary D, p → 1).** The duality π_p = π_q is a theorem of the literature (Keller–Vakil, Amer. Math. Monthly 116 (2009); the paper cites it), not the numerical check of D12 §1.2.
- Wording: the signs have period four for all sufficiently large k (numerically for k = 7 and every 9 ≤ k ≤ 150; not at k = 3, 4, 5, 6, 8). The O-constants are not computed. For even k the bracket error is in fact O(log k/k). U_k is "of order" (k−2)!/(2π)^{k−1}, not "≈". The optimal-truncation error is e^{−πp} up to a power of p.
- Referee's independent data: N_k for k ≤ 150 in certified ball arithmetic by a route that does not use Theorem A; residuals of the two-term law bounded (≤ 2.62 odd, ≤ 6.65 even, in the k²/log²k scaling).

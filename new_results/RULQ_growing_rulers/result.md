# RULQ — rulers that grow with p; Borel summability and optimal truncation of π_p

Labels: **PROVED** (complete argument below, or a complete adaptation of a refereed proof with every change written out),
**DERIVED** (mechanism clear, one named step not justified), **NUMERICAL** (computed, not proved).
Notation: h = 1/p, ε = h, ψ = Γ'/Γ, γ = Euler's constant, y₀ = iπ, Log = principal logarithm (Log y₀ = log π + iπ/2).
L_q(p) = ℓ^q-length of the ℓ^p unit circle; L_q(p)/8 = 1 + Σ_k N_k h^{k+1} (formal); for q = p this is π_p/4 (the N_k of
Theorem `thm:pip`). Every number quoted was printed by a script in this folder (logs listed in Sec. 8).

---

## 0. Verdicts in one paragraph

**A.** Bounded rulers that depend analytically on 1/p give a convergent series (PROVED). The dual ruler is exact:
L_{p/(p−1)}(p) = (4/p)B(1/p,1/p) = 8Γ(1+1/p)²/Γ(1+2/p) = 2·area of the p-disc; its series converges exactly for p > 1
and has 2ζ(3) at order p^{−3} (PROVED). Every ruler that grows linearly, q = λp + c + O(1/p) analytic in 1/p (λ > 0;
this contains q = λp and q = p + c), gives a series that diverges for every p, with
N_k = 2Re{e^{iπ(1+λ−c)}λ^{−k}[(k−2)! y₀^{1−k} + λ(k−3)! y₀^{2−k}(ψ(k−2)+γ−Log y₀ − b₁y₀ − c/λ)]} + O((k−2)!(πλ)^{−k}k^{−2}log²k)
(b₁ defined in Thm A3), and Borel singularities at t = ±iπλ (PROVED). Checked to k = 90 for λ = 2, 1/2, 3/2 and q = p+1, p+½.
Faster growth (q = p²) diverges too, with Gevrey order 1/2 (DERIVED + NUMERICAL). So, within rulers analytic in 1/p:
**the series converges iff q stays bounded** (proved for bounded and linearly growing q; derived for faster growth).
**B.** The series of π_p is Borel summable along the positive p-direction, its Borel transform is analytic in the strip
|Im t| < π with radius exactly π, and its Borel sum is π_p/4 − 1 for every p > (1+√5)/2 (PROVED). An explicit rigorous
truncation bound is proved (Theorem B2); it is informative for p ≥ 5 (above the truth by a factor 578 at p = 5, 19 at p = 10, 15 at p = 20). The true size of the
error at the smallest term is √2·p^{−5/2}e^{−πp}·(1 + O(log p/p)) (DERIVED; numerically 0.62–0.86 times that for 4 ≤ p ≤ 14, 0.89 at p = 20).
**C.** The leading beyond-all-orders term is R_K(p) ≈ √2 cos(π(K−1)/2 + π/4)·Γ(K−1)/(p²(πp)^{K−1}) (sign pattern + − − + for
K ≡ 1,2,3,0 mod 4); a refined model (the exact Borel remainder of the leading singular part) reproduces R_K to 1–4 % at
p = 8, 10, 12 (DERIVED + NUMERICAL).

---

## 1. Status table

| # | Claim | Status | Where | Numerical check |
|---|---|---|---|---|
| A0 | Exact forms of L_q/8 (u-, v-, s-forms; Laplace form for general q) | PROVED | §2 | two quadratures agree to working precision (30 dps: 2·10⁻³¹); λ=1 coefficients = PIPDIV to 4·10⁻⁹⁸ |
| A1 | Dual ruler: L_{p′}(p) = (4/p)B(1/p,1/p) = 8Γ(1+h)²/Γ(1+2h); series = 8exp(Σ_{k≥2}(−1)^k(2−2^k)ζ(k)h^k/k), radius exactly 1 | PROVED | §3.1 | quadrature vs Γ form 8·10⁻⁶² (p = 1.5, 2, 3, 7); D₂ = −ζ(2), D₃ = 2ζ(3), D₄ = −π⁴/40 to 25 digits |
| A2 | Every bounded ruler q = Q(1/p), Q analytic, Q(0) ≥ 1: series converges (radius ≥ explicit domain) | PROVED | §3.2 | fixed q: PIPDIV/REF_PIP; dual: A1 |
| A3 | q = pQ(1/p), Q analytic, Q(0) = λ > 0: asymptotic series exists, N_k two-term law, |N_k|^{1/k} ~ k/(eπλ), radius 0, Borel radius exactly πλ | PROVED (adaptation of PIPDIV, all changes listed) | §3.3 | k ≤ 90, λ = 2, ½, 3/2, q = p+1, p+½ (Table A) |
| A3′ | Same family: Borel summable along p > 0 (Borel transform analytic in |Im t| < πλ) | PROVED (same proof as B1) | §4.1 Remark | — |
| A4 | q = p^d, d ≥ 2: divergent, Gevrey order 1/d | DERIVED | §3.5 | q = p²: |N_k|^{2/k}/(k/2eπ) ∈ [1.06, 1.19], k = 30..90 |
| A5 | Optimal-truncation error for q = λp is ≈ smallest term, ≪ e^{−πλp} | NUMERICAL | §3.6 | 24 cases (Table A′) |
| B1 | π_p: Borel summable along p>0; Borel transform 𝔅(t) = Σ N_k t^k/k! analytic in |Im t| < π, radius exactly π; π_p/4 − 1 = ∫₀^∞ e^{−pt}𝔅(t)dt for p > φ = (1+√5)/2 | PROVED | §4.1 | — |
| B2 | Explicit rigorous bound |R_K(p)| ≤ R0+R1+R2 (closed formulas + one-dimensional integrals) | PROVED | §4.2 | Table B: bound/|R| = 1.1·10⁵, 4.8·10⁴, 578, 19.4, 14.9 at p = 2, 3, 5, 10, 20 |
| B3 | True size at the smallest term: |R_K| = √2 p^{−5/2}e^{−πp}e^{ν²/(2πp)}(1+O(log p/p)), ν = K−1−πp | DERIVED | §4.3 | Table C2: ratio 0.62–0.86 (even K), parity correction O(log p/p) |
| C1 | Leading beyond-all-orders term, sign pattern of period 4 in K | DERIVED | §5 | Table C |
| C2 | R_K = R_K^{sing}(1 + O(1/p)), R^{sing} = exact remainder of the l = 0 singular model | DERIVED | §5 | R/R^{sing} = 1.009–1.040 for p = 8, 10, 12 |

---

## 2. Exact forms (PROVED)

D12c §1 (refereed in REF_PIP_referee) gives, for p > 1, q ≥ 1, with r := q(1−h) and η := 1/q,

  L_q(p)/8 = 2^{−h} + Ê_q(h),  Ê_q(h) = h∫₀¹ u^{h−1}(1+u)^{−1−h}[(1+u^r)^η − 1]du.          (2.1)

*v-form.* u = e^{−v}: Ê_q = h∫₀^∞ e^{−hv}(1+e^{−v})^{−1−h}[(1+e^{−rv})^η − 1]dv.            (2.2)
*s-form.* u = e^{−κs}, κ = 1/r: Ê_q = hκ∫₀^∞ G(κs,h)[(1+e^{−s})^η − 1]ds, G(y,h) = e^{−hy}(1+e^{−y})^{−1−h}.   (2.3)

**Laplace form.** For 0 < η ≤ 1 and v > 0, (1+e^{−rv})^η − 1 = Σ_{n≥1} C(η,n)e^{−nrv}, C(η,n) = binomial(η,n), and
Σ_n|C(η,n)| = 1 (signs (−1)^{n−1}). Since e^{−hv}(1+e^{−v})^{−1−h} = e^{v}(1+e^{v})^{−1−h}, Fubini (dominating function
h e^{−hv}(1+e^{−v})^{−1−h}, integrable) gives

  Ê_q = h Σ_{n≥1} C(η,n) ∫₀^∞ e^{−nrv} e^{v}(1+e^v)^{−1−h} dv.                              (2.4)

Now let q = Q(h)/h with Q analytic near 0, Q(0) = λ > 0 (q = λp + c + …, c = Q′(0)). Then
r = λ/h + β(h) with β(h) := (Q(h)(1−h) − λ)/h analytic, b₀ := β(0) = c − λ, b₁ := β′(0) = Q″(0)/2 − c, and

  Ê_q = h Σ_{n≥1} C(η,n) ∫₀^∞ e^{−nλv/h} H_n(v,h) dv,  H_n(v,h) = e^{(1−nβ(h))v}(1+e^v)^{−1−h},  η = h/Q(h).   (2.5)

For λ = 1, Q ≡ 1 (q = p): β ≡ −1 and H_n = e^{(n+1)v}(1+e^v)^{−1−h}, which is PIPDIV Lemma 1. For q = λp: β ≡ −λ.
For q = p + c: β = c − 1 − ch. The numerical coefficients below use (2.3) expanded term by term (script ru1; J_{i,j} =
∫₀^∞ s^j log^i(1+e^{−s})ds); the proofs use (2.5). For λ = 1 the two routes agree with PIPDIV's N_k to 4·10⁻⁹⁸.

---

## 3. Task A — rulers that grow (or do not grow) with p

### 3.1 Theorem A1 (dual ruler) — PROVED

For p > 1 and p′ = p/(p−1):
  L_{p′}(p) = (4/p) B(1/p, 1/p) = 8 Γ(1+1/p)² / Γ(1+2/p) = 2·Area{|x|^p + |y|^p ≤ 1}.
Hence L_{p′}(p)/8 = exp(Σ_{k≥2} (−1)^k (2−2^k) ζ(k) p^{−k}/k) = 1 − ζ(2)p^{−2} + 2ζ(3)p^{−3} − (π⁴/40)p^{−4} + …,
and this series in 1/p converges exactly for p > 1 (radius 1 in h; simple pole at h = −1).

*Proof.* In L_q(p) = (4/p)∫₀¹(t^{q/p−q} + (1−t)^{q/p−q})^{1/q}dt put q = p′: q/p − q = q(1/p − 1) = −1 and 1/q = 1 − 1/p.
So the integrand is (1/t + 1/(1−t))^{1−1/p} = (t(1−t))^{1/p−1}, and the integral is B(1/p,1/p). Then
(4/p)Γ(h)²/Γ(2h) = 8Γ(1+h)²/Γ(1+2h) (use Γ(1+x) = xΓ(x)); the area of the p-disc is 4Γ(1+h)²/Γ(1+2h).
For |h| < 1, log Γ(1+h) = −γh + Σ_{k≥2}ζ(k)(−h)^k/k, so 2logΓ(1+h) − logΓ(1+2h) = Σ_{k≥2}(−1)^k(2−2^k)ζ(k)h^k/k (the γ
terms cancel); this converges for |h| < 1/2 and gives the displayed series there. The function Γ(1+h)²/Γ(1+2h) is
meromorphic; its only possible poles are at h = −1, −2, … (Γ(1+h)² has double poles there, 1/Γ(1+2h) is entire with a simple
zero at h = −1). At h = −1 the order is 2 − 1 = 1, so there is a simple pole; it is the singularity nearest to 0, and the
Taylor series at 0 has radius exactly 1. □
(Remark: "length of the unit circle of a norm, measured in the dual norm turned by 90°, equals twice the area" is a classical
fact of Minkowski geometry; the ℓ^{p′} ball is invariant under the 90° turn. We do not need it.)
*Numerical:* quadrature of (2.2) vs 8Γ(1+h)²/Γ(1+2h) at p = 1.5, 2, 3, 7: 7.8·10⁻⁶² (60 dps). 80-term partial sums:
error 0 (p=10, 60 dps), 5.5·10⁻²⁵ (p=2), 1.6·10⁻⁸ (p=1.25, slow, as radius 1 predicts); |D_k| → 2 (pole residue). [ru2_log]

### 3.2 Theorem A2 (bounded analytic rulers converge) — PROVED

Let Q be analytic on |h| < ρ_Q, real on real h, with Q(0) ≥ 1, and let q(p) = Q(1/p). Let U be the connected component
containing 0 of {|h| < ρ_Q : Q(h) ≠ 0, Re[Q(h)(1−h)] > 0, Re h + Re[Q(h)(1−h)] > 0}. Then F(h) := L_{q(1/h)}(1/h)/8,
defined for small h > 0, extends analytically to U. Hence its 1/p series converges, with radius at least dist(0, ∂U) > 0.
For fixed q ≥ 1, U ⊃ {|h| < 1} (convergence for p > 1, as in prop:ruler). For the dual ruler U ⊃ {|h| < 1}.

*Proof.* Use (2.1): F(h) = 2^{−h} + h∫₀¹ u^{h−1}(1+u)^{−1−h}Ψ(u,h)du, Ψ = (1+z)^{1/Q(h)} − 1, z = u^{r(h)} = e^{r(h)log u},
r(h) = Q(h)(1−h), the power defined with the principal Log(1+z). For h ∈ U and u ∈ (0,1): |z| = u^{Re r} < 1, so
Re(1+z) > 0 and Ψ is analytic in h. For h in a compact K ⊂ U: the set {u^{r(h)} : u ∈ [0,1], h ∈ K} is compact and does not
contain −1 (|z| < 1 for u < 1, z = 1 at u = 1), so |Log(1+z)| ≤ C_K|z| and, since |1/Q| is bounded on K,
|Ψ| ≤ |Log(1+z)/Q|·e^{|Log(1+z)/Q|} ≤ C′_K u^{Re r(h)}. The integrand is then bounded by
C″_K u^{Re h + Re r(h) − 1}, and Re h + Re r(h) ≥ σ_K > 0 on K, so it is dominated by an integrable function uniformly on K.
Holomorphy of each integrand in h, local uniform domination and Morera's theorem (Fubini on triangles) make F analytic on U.
F equals L/8 for small real h > 0, so the Taylor series of F at 0 is the 1/p expansion, and it converges in the largest disc
inside U. For fixed q: r = q(1−h), Re r > 0 ⇔ Re h < 1, and Re h + q(1 − Re h) > 0 holds for Re h < 1 (q ≥ 1); so U ⊃ {Re h < 1}.
Dual: r ≡ 1, Q = 1/(1−h) ≠ 0, U ⊃ {Re h > −1, h ≠ 1} ⊃ {|h| < 1}. □

### 3.3 Theorem A3 (linear growth: divergence and growth law) — PROVED

Let Q be analytic near 0, real on reals, Q(0) = λ > 0, and q(p) = pQ(1/p) = λp + c + O(1/p) (c = Q′(0)). With b₀, b₁ as in
§2 (b₀ = c − λ, b₁ = Q″(0)/2 − c), y₀ = iπ:

(i) L_{q(p)}(p)/8 has an asymptotic expansion 1 + Σ_{k≥0}N_k h^{k+1} to all orders, and
  N_k = (−log 2)^{k+1}/(k+1)! + T_k + U_k,  T_k = Σ_{m=0}^{k−2} m! λ^{−m−1}[h^{k−2−m}](a_m(h)/Q(h)),
  a_m(h) = [v^m]A(v,h), A(v,h) := H₁(v,h) = e^{(1−β(h))v}(1+e^v)^{−1−h},  |U_k| ≤ C(2^k + k²(k−2)!(1.9πλ)^{−k}).
(ii) As k → ∞,
  N_k = 2Re{ e^{(1−b₀)y₀} λ^{−k} [ (k−2)! y₀^{1−k} + λ(k−3)! y₀^{2−k}(ψ(k−2) + γ − Log y₀ − b₁y₀ − c/λ) ] }
        + O((k−2)! π^{1−k} λ^{−k} log²k / k²).
  In real form, with θ = π(1+λ−c) and X_k = ψ(k−2)+γ−log π:
  N_k = (2(k−2)!/(π^{k−1}λ^k)) cos(θ − π(k−1)/2)
      + (2(k−3)!/(π^{k−2}λ^{k−1})) [ (X_k − c/λ) cos(θ − π(k−2)/2) + π(½ + b₁) sin(θ − π(k−2)/2) ] + O(…)
  (from (X_k − iπ/2 − iπb₁ − c/λ) and Re[e^{iα}(u − iw)] = u cos α + w sin α). For linear Q, b₁ = −c and the bracket is
  (X_k − c/λ)cos(·) − π(c − ½)sin(·). Check (λ = 1, c = 0, odd k): it gives the PIPDIV correction factor 1 + π²/(2(k−2)).
(iii) The leading term never vanishes for two consecutive k, so |N_k|^{1/k} = (k/(eπλ))(1+o(1)) along a subsequence and
  limsup |N_k|^{1/k}/k = 1/(eπλ). The series diverges for every p (radius 0).
(iv) The Borel transform Σ_k N_k t^k/k! has radius of convergence exactly πλ; its singularities on |t| = πλ are at ±iπλ.

Special cases (c = 0 for q = λp; λ = 1 for q = p + c):
- q = λp: N_k = −(2(k−2)!/(π^{k−1}λ^k))·[cos(π(λ − (k−1)/2)) + O(log k/k)]; for λ ∈ ℤ the odd k carry the leading term, for
  λ ∈ ℤ+½ the even k do; at the other parity the first term vanishes and N_k is of size (k−3)!log k/(π^{k−2}λ^{k−1}).
- q = p + c: N_k = (2(k−2)!/π^{k−1})·[cos(π(c + (k−1)/2)) + O(log k/k)]. For c = 1 the sign pattern of the odd k flips
  relative to π_p; for c = ½ the even k lead and the O(log k/k) term vanishes identically at even k (checked: ν₉₀ = −0.99476
  vs −1).
- λ = 1, c = 0 is PIPDIV Theorem C.

*Proof.* It is PIPDIV §2 with the following changes; everything else is word for word.
1. *Laplace form* (PIPDIV Lemma 1 → (2.5)). Valid for 0 < h < h₀ (so that 0 < η = h/Q(h) ≤ 1 and r > 0); proof as in §2.
2. *Watson step, uniform in n* (Lemma 2). Put t = nλv: nλ∫₀^∞e^{−nλv/h}H_n dv = ∫₀^∞e^{−t/h}H_n(t/(nλ),h)dt. The symbol
   H_n(t/(nλ),h) = e^{t/(nλ)}e^{−β(h)t/λ}(1+e^{t/(nλ)})^{−1−h} is analytic in |t| < πλ for every n ≥ 1 (and in |t| < 2πλ for
   n ≥ 2), and uniformly bounded for |t| ≤ 0.95πλ (n ≥ 1), |t| ≤ 1.9πλ (n ≥ 2), |h| ≤ h₁ (β analytic and bounded there). Cauchy
   gives |b_{n,m}| ≤ M₀(0.95πλ)^{−m}, resp. M₁(1.9πλ)^{−m}. Tail: for real 0 ≤ h ≤ h₁, H_n(v,h) ≤ e^{−(nβ(h)+h)v}, so
   H_n(t/(nλ),h) ≤ e^{B t} with B = sup|β|/λ; the integral over [T₀,∞) is ≤ 2he^{−T₀/(2h)} for h ≤ 1/(2B). Binomial: for
   |η| ≤ ½, |C(η,n)| ≤ e^{1/2}|η|n^{−1/2}, and |η| ≤ 2|h|/λ for small h. The rest of Lemma 2 is unchanged and gives (i); the
   n = 1 part is T_k (C(η,1) = h/Q(h), b_{1,m} = λ^{−m}a_m). U_k is bounded as in PIPDIV with 6 replaced by 1.9πλ, and is
   O(k²(π/1.9π)^k) = o(k^{−L}) relative to (k−2)!(πλ)^{−k}.
3. *Local structure* (Lemma 3). 1+e^v = (y₀ − v)φ(v − y₀) (φ(u) = (e^u−1)/u) since e^{y₀} = −1, so on V_± (and on the larger
   connected set W of the referee's repair R1)
   A(v,h) = e^{(1−β(h))y₀}(y₀ − v)^{−1−h}Φ(v − y₀,h),  Φ(u,h) = e^{(1−β(h))u}φ(u)^{−1−h}.
   The branch argument is the same: the factor e^{(1−β(h))v} is entire and single-valued.
4. *Hankel bound* (Lemma 4). At h = 0, A(·,0) = e^{(1−b₀)v}/(1+e^v) has simple poles at ±iπ with residue −e^{(1−b₀)y₀}, which is
   also the residue of the l = 0 term e^{(1−b₀)y₀}(y₀−v)^{−1}; so R_L(·,0) is analytic in |v| < 3π. For |h| ≤ h₁ the local
   bound uses max |Φ| on |u| = π, which is finite because e^{(1−β(h))u} is bounded there. No other change.
5. *Proposition 5 and Theorem C.* T_k has the extra weights λ^{−m−1} and the analytic factor 1/Q(h). In Prop. 5 each term
   m!λ^{−m−1}[h^i][y^m](R_L/Q) is bounded by e C_L m!(log k)^i m^{−L}(πλ)^{−m}·sup_{|h|≤1/log k}|1/Q|, so the same summation
   gives T_k = T_k^{(L)} + O((k−2)!(πλ)^{−k}k^{−L−1}log k) with
   T_k^{(L)} = 2Re Σ_{l<L}Σ_m λ^{−m−1}[h^{k−2−m}]((−1)^lφ_l(h)e^{(1−β(h))y₀}y₀^{l−1−m}e^{−hLog y₀}(1+h−l)_m/Q(h)).
   For L = 1, with F_M(h) = e^{(1−β(h))y₀ − hLog y₀}Γ(M+h)/(Γ(1+h)Q(h)), M = m+1: [h⁰]F_M = e^{(1−b₀)y₀}Γ(M)/λ,
   [h¹]F_M = e^{(1−b₀)y₀}(Γ(M)/λ)(ψ(M) + γ − Log y₀ − b₁y₀ − c/λ) (because (log Q)′(0) = c/λ and (1−β)′(0) = −b₁).
   The terms i = 0, 1 (i = k−2−m) give the displayed formula; i ≥ 2 are bounded by the Gamma bound of PIPDIV §2.4 (the extra
   factors are bounded on |h| ≤ 1/log k) and give O((k−2)!(πλ)^{−k}λ·log²k/k²); Theorem B with L = 1 adds O(k^{−2}log k)
   relative. Real form: insert y₀^{1−k} = π^{1−k}e^{−iπ(k−1)/2} and take real parts.
6. *(iii)* θ − π(k−1)/2 and θ − πk/2 differ by π/2, so max(|cos|) over two consecutive k is ≥ 1/√2; Stirling.
7. *(iv)* |N_k|/k! ≍ k^{−2}(πλ)^{−k} along that subsequence and ≲ k^{−2}(πλ)^{−k}log k always, so the radius is πλ. The
   Borel transform is analytic in |Im t| < πλ (Remark after B1), so the only points of |t| = πλ where it can be singular are
   ±iπλ, and by conjugate symmetry both are singular. □

**Interpretation.** The Borel singularities sit at t = ±iπλ in the variable conjugate to p, i.e. at ±iπ in the variable
conjugate to the ruler exponent q ≈ λp: the divergence is a "1/q" effect, with the same pair of zeros of 1+e^v at ±iπ as for π_p.
The phase e^{iπ(1+λ−c)} = −e^{iπ(λ−c)} only decides which parity of k carries the factorial growth.

### 3.4 The general picture (which rulers converge)

Within rulers analytic in 1/p at p = ∞ (q = p^d Q(1/p), Q analytic, Q(0) > 0, d ∈ {0,1,2,…}):
d = 0 (bounded): convergent (A2, PROVED); d = 1: divergent, Gevrey 1, Borel radius πλ (A3, PROVED); d ≥ 2: divergent,
Gevrey order 1/d (A4, DERIVED). So **the series converges iff q stays bounded** (PROVED when q is bounded or grows
linearly; for faster polynomial growth the divergence is DERIVED). Rulers that are not analytic in 1/p
(q = √p, q = log p) do not have a power series in 1/p at all (the expansion contains fractional powers or logs), so the question
does not arise in the same form; the same mechanism (kernel e^{−rv} with r → ∞) still produces Gevrey divergence in the
appropriate variable (DERIVED, not checked).

### 3.5 A4: q = p^d, d ≥ 2 — DERIVED + NUMERICAL

For q = Q/h^d, the n = 1 term of (2.4) is hη∫₀^∞e^{−rv}e^{v}(1+e^v)^{−1−h}dv with r = Q(1−h)/h^d. With w = rv it is
(hη/r)Σ_m m! r^{−m}[v^m](e^{v}(1+e^v)^{−1−h}) (Watson), and [v^m] ~ π^{−m} (poles at ±iπ). Since r^{−1} = O(h^d), the terms grow
like m!(h^d/π)^m: Gevrey order 1/d, Borel singularities at p^d t = ±iπ/Q(0). Not proved: the re-expansion of
(1−h)^{−m−1} and of the h-dependent factors mixes many m into one power of h, and a lower bound on |N_k| would need the
PIPDIV machinery in two scales. *Numerical* (q = p², k ≤ 90): |N_k|^{2/k}/(k/(2eπ)) = 1.06, 1.19, 1.19, 1.16, 1.07, 1.14,
1.13 at k = 30, 40, …, 90: consistent with |N_k| = (k/(2eπ))^{k/2}e^{O(k)}, i.e. Gevrey order exactly 1/2. [ru5_log]

### 3.6 Numerical checks for A

Coefficients N_k (k ≤ 90) were computed from (2.3) with two independent J-tables (100 dps, 192-node panels, S = 700;
80 dps, 96-node panels, S = 640): runs agree to ≤ 4.7·10⁻⁷⁷ relative for λ = 2, ½, q = p+1, q = p. The λ = 1 coefficients equal
PIPDIV's to 4·10⁻⁹⁸. [ru1_crosscheck_log]

**Table A** (ν_k = N_k/(2(k−2)!π^{1−k}λ^{−k}); F = two-term law; lead = cos(θ − π(k−1)/2); r_k = (ν_k − F_k)k²/log²k,
which the theorem says is bounded; last column → 1 for the parity that carries the leading term):

| ruler | k | ν_k | F_k | lead | r_k | |N_k|^{1/k}eπλ/k |
|---|---|---|---|---|---|---|
| q = 2p | 89 | −1.08769 | −1.11344 | −1 | 10.1 | 0.957 |
| q = 2p | 90 | −0.31147 | −0.27875 | 0 | −13.1 | 0.944 |
| q = p/2 | 89 | 0.07233 | 0.07028 | 0 | 0.81 | 0.929 |
| q = p/2 | 90 | −1.02686 | −1.02804 | −1 | 0.47 | 0.957 |
| q = 3p/2 | 89 | −0.22953 | −0.21085 | 0 | −7.3 | 0.941 |
| q = 3p/2 | 90 | 1.07060 | 1.08412 | 1 | −5.4 | 0.957 |
| q = p+1 | 89 | −0.94072 | −0.94328 | −1 | 1.01 | 0.956 |
| q = p+1 | 90 | −0.09793 | −0.10368 | 0 | 2.30 | 0.932 |
| q = p+½ | 89 | 0.12248 | 0.12251 | 0 | −0.01 | 0.934 |
| q = p+½ | 90 | −0.99476 | −1.00000 | −1 | 2.10 | 0.957 |
| q = p | 89 | 1.05108 | 1.05672 | 1 | −2.2 | 0.957 |
| q = p | 90 | 0.14746 | 0.13938 | 0 | 3.2 | 0.937 |

More rows (k = 30, 31, 50, 51, 70, 71, 88) are in ru5_log.txt: r_k stays in [−25, 13] for all families. A sharper test of the
singular structure: with S_k the full l = 0 singular part (all orders in h), (N_k − S_k)k²/scale at k = 73 and 89 is
−57.5, −56.1 (λ = 2); −4.86, −4.88 (q = p+1); 0.51, 0.43 (λ = ½, the suppressed parity): N_k − S_k = O(scale/k²), as
Lemma 4 with L = 1 predicts.
T_k/N_k − 1 = 1.7·10⁻²⁷ (λ = 2, k = 89): the n = 1 term carries everything, as in (i). [ru1_log_*]

**Table A′** (direct check that the series is the expansion of the true function): |L/8 − S_{K*}| with S_{K*} the sum up to
the smallest term, L/8 by quadrature of (2.2) (two quadratures agree to working precision):

| ruler | p = 4 | p = 8 | p = 12 | p = 16 | e^{−πλp} at p = 16 |
|---|---|---|---|---|---|
| q = 2p | 5.5e−14 | 2.5e−25 | 1.3e−36 | 1.6e−47 | 2.2e−44 |
| q = p/2 | 6.0e−4 | 7.1e−8 | 4.7e−11 | 4.3e−14 | 1.2e−11 |
| q = 3p/2 | 6.1e−11 | 1.3e−19 | 3.3e−28 | 1.1e−36 | 1.8e−33 |
| q = p+1 | 7.6e−8 | 5.9e−14 | 8.6e−20 | 1.6e−25 | 1.5e−22 |
| q = p | 1.2e−7 | 7.6e−14 | 9.9e−20 | 1.7e−25 | 1.5e−22 |

The error is of the size of the smallest term (ratio 0.1–3 except q = p+½ at p = 4, where an accidentally small early term
stops the sum too soon), and the smallest term is ≈ e^{−πλp} times a negative power of p. [ru2_log]

---

## 4. Task B — Borel summability and optimal truncation for π_p

Throughout q = p, π_p/4 = 2^{−h} + Ê(h), Ê = Σ_{n≥1}T_n, T_n(h) = hC(h,n)∫₀^∞e^{−ny/h}H_n(y,h)dy,
H_n = e^{(n+1)y}(1+e^y)^{−1−h} (PIPDIV Lemma 1, 0 < h < 1).

### 4.1 Theorem B1 (Borel summability) — PROVED

Let 𝔅(t) = Σ_{k≥0} N_k t^k/k!. Then (a) 𝔅 has radius of convergence exactly π; (b) 𝔅 extends analytically to the strip
S = {|Im t| < π}, and for every ε > 0 and every closed sub-strip, |𝔅(t)| ≤ C(1+|t|)e^{(φ+ε)|t|} there (φ = (1+√5)/2);
(c) for every p > φ, π_p/4 − 1 = ∫₀^∞ e^{−pt}𝔅(t)dt (absolutely convergent). So the series Σ N_k p^{−k−1} is Borel summable
in the direction of positive p, and its Borel sum is π_p/4 − 1. The only singularities of 𝔅 on |t| = π are ±iπ.

*Proof.* Step 1 (parametric Laplace form). Put t = ny in T_n: T_n = (hC(h,n)/n)∫₀^∞e^{−t/h}H_n(t/n,h)dt. Define
  B(t,h) := Σ_{n≥1} (hC(h,n)/n) H_n(t/n,h).
For real 0 < h < 1 and t ≥ 0, 0 ≤ H_n(t/n,h) ≤ e^{(n+1)t/n}e^{−(1+h)t/n} ≤ e^t and Σ|C(h,n)| = 1, so Fubini gives
Ê(h) = ∫₀^∞e^{−t/h}B(t,h)dt for 0 < h < 1.
Step 2 (analyticity and bounds). Fix 0 < ρ < 1 and 0 < δ ≤ π/2; let S_δ = {|Im t| ≤ π − δ}. For |η| ≤ ρ, n ≥ 1, t ∈ S_δ and
v = t/n (so |Im v| ≤ π − δ): (a) |1+e^v|² = 1 + 2e^x cos y + e^{2x} ≥ (1 − cos δ)(1+e^{2x}) (x = Re v, y = Im v; because
cos y ≥ −cos δ and cos δ(1−e^x)² ≥ 0), so |1+e^v| ≥ σ√(1+e^{2x}), σ = √2 sin(δ/2) ≤ 1; (b) 1+e^v ∉ (−∞,0], so the principal
power is analytic and |(1+e^v)^{−1−η}| ≤ |1+e^v|^{−1−Re η}e^{π|Im η|} ≤ σ^{−1−ρ}(1+e^{2x})^{−(1−ρ)/2}e^{πρ}. Hence
|H_n(t/n,η)| ≤ e^{πρ}σ^{−1−ρ}e^{(n+1)x − (1−ρ)x⁺} ≤ e^{πρ}σ^{−1−ρ}e^{(1+ρ)(Re t)⁺}. And |C(η,n)| = (|η|/n)Π_{i<n}|1−η/i| ≤
(ρ/n)e^{ρH_{n−1}} ≤ ρe^ρ n^{ρ−1}. So the series for B converges uniformly on S_δ × {|η| ≤ ρ}:
  |B(t,η)| ≤ C_{ρ,δ} e^{(1+ρ)(Re t)⁺},  C_{ρ,δ} = ρ²e^{ρ(1+π)}σ^{−1−ρ}ζ(2−ρ).
Each term is analytic, so B is analytic on S × {|η| < 1}.
Step 3 (Lemma P, Borel transform of a parametric Laplace integral). Write B(t,η) = Σ_{j≥0}B_j(t)η^j; by Cauchy,
|B_j(t)| ≤ C_{ρ,δ}ρ^{−j}e^{(1+ρ)(Re t)⁺} on S_δ. Let (I^j f)(t) = ∫₀^t (t−τ)^{j−1}/(j−1)! f(τ)dτ (straight segment; I⁰f = f) and
B̂ := Σ_j I^jB_j. For t ∈ S_δ the segment [0,t] lies in S_δ, and parametrising it by arc length u ∈ [0,|t|],
  |I^jB_j(t)| ≤ Cρ^{−j}∫₀^{|t|}(|t|−u)^{j−1}/(j−1)!·e^{(1+ρ)u}du,  so  Σ_{j≥1}|I^jB_j(t)| ≤ Cρ^{−1}∫₀^{|t|}e^{(|t|−u)/ρ}e^{(1+ρ)u}du;
the series converges locally uniformly, B̂ is analytic on S, and |B̂(t)| ≤ C(1 + |t|/ρ)e^{max(1+ρ,1/ρ)|t|} on S_δ.
For real h with 0 < h < ρ and 1/h > 1+ρ: ∫₀^∞e^{−t/h}I^jB_j dt = h^j∫₀^∞e^{−t/h}B_j dt (Laplace transform of a convolution with
t^{j−1}/(j−1)!; absolutely convergent), and Σ_j h^j∫e^{−t/h}|B_j|dt ≤ Σ_j(h/ρ)^j C/(1/h − 1 − ρ) < ∞, so
  ∫₀^∞e^{−t/h}B̂(t)dt = Σ_j h^j∫₀^∞e^{−t/h}B_j(t)dt = ∫₀^∞e^{−t/h}B(t,h)dt = Ê(h).
B̂ does not depend on ρ: its Taylor coefficients at 0 are Σ_{j+m=k}[t^m]B_j·m!/k!, built from B alone.
Step 4 (identification). B̂ is analytic at 0 and exponentially bounded on [0,∞); Watson's lemma gives
∫e^{−t/h}B̂ ~ Σ_k k![t^k]B̂·h^{k+1}. By uniqueness of asymptotic expansions and PIPDIV Theorem A,
k![t^k]B̂ = N_k − (−log2)^{k+1}/(k+1)!. The convergent part: β₀(t) := Σ_k (−log2)^{k+1}t^k/(k!(k+1)!) is entire,
|β₀(t)| ≤ log2·e^{2√(log2·|t|)}, and ∫₀^∞e^{−t/h}β₀ = Σ(−h log2)^{k+1}/(k+1)! = 2^{−h} − 1 (termwise, dominated). So
𝔅 = β₀ + B̂: this gives (b) (choose ρ ∈ (1/φ − ε′, 1/φ)), and (c): for p > φ choose ρ with 1/p < ρ < p − 1 (possible iff
p² − p − 1 > 0); then h = 1/p < ρ, 1/h = p > 1+ρ, and p > max(1+ρ,1/ρ), so the Laplace integral converges absolutely and
equals 2^{−h} − 1 + Ê(h) = π_p/4 − 1.
Step 5 (a). By PIPDIV Theorem C, |N_k|/k! ≍ k^{−2}π^{−k} for odd k and ≲ k^{−2}π^{−k}log k, so the radius is π. 𝔅 is analytic
in |t| < π and in S; the points of |t| = π outside S are ±iπ; at least one is singular, and 𝔅(t̄) = conj 𝔅(t), so both are. □

*Remark (A3′, growing rulers).* For q = pQ(1/p) as in A3, the same proof applies to (2.5): B(t,h) = Σ_n (hC(η,n)/(nλ))
H_n(t/(nλ),h) is analytic for |Im t| < πλ, |h| small, and exponentially bounded (the factor e^{−β(h)t/λ} adds a fixed
exponential rate). So those series are Borel summable along p > 0 too (for p larger than an explicit constant), with Borel
transform analytic in |Im t| < πλ and radius exactly πλ.

*Remark (what Borel summability does not say).* The singularities ±iπ are off the positive axis, so there is no Stokes
ambiguity and no exponentially small term is "missing" from the Borel sum. The e^{−πp} that appears below is only the size of
the optimally truncated remainder.

### 4.2 Theorem B2 (explicit rigorous truncation bound) — PROVED

Let R_K(p) := π_p/4 − 1 − Σ_{k<K}N_k p^{−k−1}, h = 1/p, K ≥ 5.

**Lemma R (remainder of a Laplace transform).** Let 0 < r < π′, f analytic on a neighbourhood of S_r = {|Im s| ≤ r}, real on
R, |f(s)| ≤ Me^{a(Re s)⁺} on S_r, and |f(x ± ir)| ≤ m(x) (x ∈ R). For N ≥ 2 and 0 < h < 1/a let
Rem_N[f](h) = ∫₀^∞e^{−y/h}f(y)dy − Σ_{m<N}m! f_m h^{m+1} (f_m Taylor coefficients at 0). Then
  |Rem_N[f](h)| ≤ (1/(πr))∫_ℝ m(x)(x²+r²)^{−N/2} h^{N+1}Γ(N+1,(x−r)⁺/h) dx + (M/π)e^{ar}∫₀^∞e^{−y(1/h−a)}y^N(y+r)^{−N}dy.
*Proof.* Since ∫e^{−y/h}y^m = m!h^{m+1}, Rem_N = ∫₀^∞e^{−y/h}(f − T_Nf)(y)dy, T_N the Taylor polynomial of degree N−1. For
y ≥ 0, 1/(s−y) = Σ_{m<N}y^m/s^{m+1} + y^N/(s^N(s−y)) gives f(y) − T_Nf(y) = (y^N/2πi)∮f(s)ds/(s^N(s−y)) on the boundary of the
rectangle [−X, y+r]×[−r, r]. The left edge contributes ≤ y^N r M/(πX^{N+1}) → 0 (X → ∞). On the horizontal edges
|s − y| ≥ r, |s|^{−N} = (x²+r²)^{−N/2}, and |f| ≤ m(x) on both (f real, so |f(x−ir)| = |f(x+ir)|); they contribute
≤ (y^N/(πr))∫_{−∞}^{y+r}m(x)(x²+r²)^{−N/2}dx. On the right edge (length 2r) |s−y| ≥ r, |s| ≥ y+r, |f| ≤ Me^{a(y+r)}: it
contributes ≤ (M/π)e^{a(y+r)}y^N(y+r)^{−N}. Multiply by e^{−y/h}, integrate over y ≥ 0, and use Tonelli in the first term:
∫_{y ≥ (x−r)⁺}e^{−y/h}y^Ndy = h^{N+1}Γ(N+1,(x−r)⁺/h). □
(Corollary used for crude parts: with m(x) = Me^{ax⁺}, the bound is ≤ (M/π)(1+c_N)e^{ar}N!(h/(1−ah))^{N+1}r^{−N},
c_N = √πΓ((N−1)/2)/Γ(N/2) ≤ π.)

**Theorem B2.** Let p > 1 and choose ρ, δ, ρ₂, δ₂ with h < ρ, (1+ρ)h < 9/10, 0 < δ ≤ π/2, h < ρ₂ < 1, 0 < δ₂ < 2π. Put
r₁ = π − δ, σ₁ = √2 sin(δ/2), Λ = π + log(2/σ₁), M₁ = e^{ρΛ}/σ₁, a₁ = 1+ρ,
m₁(x) = e^{2x}e^{ρ|Log(1+e^{x+ir₁})|}/|1+e^{x+ir₁}|; r₂ = 2π − δ₂, σ₂ = √2 sin(δ₂/4), a₂ = 1+ρ₂/2,
M₂ = ρ₂²e^{ρ₂(1+π)}σ₂^{−1−ρ₂}(ζ(2−ρ₂)−1), m₂(x) = M₂e^{a₂x⁺}. Write Bnd_N(m;M,a,r) for the right side of Lemma R. Then
|R_K(p)| ≤ R0 + R1 + R2 with
  R0 = 2(h log2)^{K+1}/(K+1)!,
  R1 = h²[Σ_{j=0}^{K−4} h^jρ^{−j}Bnd_{K−2−j}(m₁;M₁,a₁,r₁) + 2M₁h(h/ρ)^{K−3}/((1−h/ρ)(1−a₁h))],
  R2 = Σ_{j=2}^{K−2} h^jρ₂^{−j}Bnd_{K−j}(m₂;M₂,a₂,r₂) + 2M₂h(h/ρ₂)^{K−1}/((1−h/ρ₂)(1−a₂h)).

*Proof.* Split R_K = R⁽⁰⁾ + R⁽¹⁾ + R⁽≥2⁾ according to π_p/4 = 2^{−h} + T₁ + Σ_{n≥2}T_n; by PIPDIV Theorem A,
N_k = (−log2)^{k+1}/(k+1)! + T_k + U_k, where T_k and U_k are the asymptotic coefficients of T₁ and of Σ_{n≥2}T_n.
(0) |Σ_{m>K}(−h log2)^m/m!| ≤ 2(h log2)^{K+1}/(K+1)! (h log2 < 1).
(1) T₁ = h²∫₀^∞e^{−y/h}A(y,h)dy, A(s,η) = g(s)e^{−ηℓ(s)}, g = e^{2s}/(1+e^s), ℓ = Log(1+e^s) (analytic on |Im s| < π). Let
A_j = [η^j]A. On S_{r₁}: |g(s)| ≤ e^{2x}/(σ₁√(1+e^{2x})) ≤ e^{x⁺}/σ₁ (Step 2 (a) of B1), and |ℓ| ≤ |log|1+e^s|| + π ≤
log(2/σ₁) + x⁺ + π = Λ + x⁺ (σ₁√(1+e^{2x}) ≤ |1+e^s| ≤ 1+e^x ≤ 2e^{x⁺}). So by Cauchy on |η| = ρ:
|A_j(s)| ≤ ρ^{−j}M₁e^{a₁x⁺} on S_{r₁}, and on the lines Im s = ±r₁, |A_j| ≤ ρ^{−j}m₁(x). For 0 < h < ρ, h < 1/a₁, dominated
convergence gives T₁ = h²Σ_j h^j∫e^{−y/h}A_j dy, and the asymptotic coefficient of h^{k+1} in T₁ is Σ_{j+m=k−2}m![y^m]A_j = T_k.
Hence R⁽¹⁾ = h²Σ_j h^jRem_{K−2−j}[A_j], with Rem_N = full Laplace transform for N ≤ 0. Lemma R for j ≤ K−4 (N ≥ 2); for
j ≥ K−3 (N ≤ 1) |Rem_N[A_j]| ≤ |∫e^{−y/h}A_j| + h|A_j(0)| ≤ 2ρ^{−j}M₁h/(1−a₁h), summed geometrically.
(≥2) Σ_{n≥2}T_n = ∫₀^∞e^{−t/h}B_{≥2}(t,h)dt with B_{≥2} = Σ_{n≥2}(hC(h,n)/n)H_n(t/n,h) (Fubini as in B1). For |η| = ρ₂, n ≥ 2,
t ∈ S_{r₂}: |Im(t/n)| ≤ π − δ₂/2, so Step 2 of B1 with δ replaced by δ₂/2 gives |H_n(t/n,η)| ≤ e^{πρ₂}σ₂^{−1−ρ₂}e^{(1+ρ₂/2)(Re t)⁺}
(n ≥ 2), and |ηC(η,n)/n| ≤ ρ₂²e^{ρ₂}n^{ρ₂−2}; summing over n ≥ 2, |B_{≥2}(t,η)| ≤ M₂e^{a₂(Re t)⁺}. B_{≥2} = O(η²). Cauchy and
the same argument give R⁽≥2⁾ = Σ_{j≥2}h^jRem_{K−j}[B_j], bounded by Lemma R with m₂ = M₂e^{a₂x⁺}. The condition (1+ρ)h < 9/10
is used only in the numerical evaluation (tail of the x-integral, below). □

*Evaluation* (ru4_bound.py): the x-integral in Lemma R is computed by quadrature on [−50r, X_R] (X_R = r + 4Nh + 40h + 10)
plus two explicit tail bounds: x < −50r: m ≤ M, Γ-factor ≤ N!h^{N+1}, giving M N!h^{N+1}(50r)^{1−N}/(N−1); x > X_R: Γ(N+1,z) ≤
e^{−θz}N!(1−θ)^{−N−1} with θ = 9/10, giving M N!h^{N+1}(1−θ)^{−N−1}X_R^{−N}e^{θr/h}e^{−(θ/h−a)X_R}/(θ/h−a). Γ(N+1,z) is
evaluated as the finite sum N!e^{−z}Σ_{i≤N}z^i/i!. The quadrature is not interval arithmetic: the bound values in Table B are
numerical evaluations of a proved formula (smooth positive integrands, 25 digits), not certified numbers.

**Asymptotic size of the bound.** The j = 0 term dominates R1 once ρ ≥ 3h (the end terms of the j-sum are of relative size ≈ (eh/ρ)^K); its
x-integral is dominated by |x| ≲ 1/√p, where m₁ ≍ (x² + δ²)^{−1/2}; with δ ≍ 1/p this gives a factor O(log p). So the bound is
O(p^{−5/2}log p·e^{−πp}) (DERIVED for the bound). Table B: bound·p^{5/2}e^{πp} = 563, 22.2, 18.8 at p = 5, 10, 20.

### 4.3 Table B (p = 2, 3, 5, 10, 20) and the true size (B3)

K_s = index of the smallest term |N_k p^{−k−1}| (k ≥ 2), truncation = Σ_{k<K_s}; K_opt = argmin_K |R_K|;
ref = √2 p^{−5/2}e^{−πp}. π_p/4 at 60 digits (v-form; checked against the w-form and π_2 = π). [ru3_log, ru4_log*]

| p | K_s | |t_{K_s}| | R_{K_s} | K_opt | R_{K_opt} | ref | bound (B2) at K_s | bound/|R_{K_s}| |
|---|---|---|---|---|---|---|---|---|
| 2 | 8 | 2.39e−4 | 3.78e−3 | 4 | 1.79e−4 | 4.67e−4 | 4.2e2 (ρ = 0.6, δ = 1.5) | 1.1e5 |
| 3 | 8 | 6.21e−6 | 8.04e−5 | 12 | 3.83e−6 | 7.32e−6 | 3.90 | 4.8e4 |
| 5 | 18 | 4.94e−9 | −2.63e−9 | 14 | −1.62e−9 | 3.81e−9 | 1.52e−6 (ρ = 0.6, δ = 0.5) | 578 |
| 10 | 34 | 6.92e−17 | −8.23e−17 | 32 | 7.86e−17 | 1.02e−16 | 1.60e−15 (ρ = 0.35, δ = 0.05) | 19.4 |
| 20 | 66 | 1.59e−31 | −3.64e−31 | — | — | 4.08e−31 | 5.41e−30 (ρ = 0.2, δ = 0.05; one admissible choice, not optimised) | 14.9 |

Reading: for p ≤ 3 the asymptotic regime has not started (K_s = 8 is far from πp + 1, the even/odd zig-zag decides the smallest
term), and B2 is not informative there. For p ≥ 5 the bound is rigorous and within the factor shown.

**B3 (true size) — DERIVED.** Model the tail by the leading large-order term 2Re[(k−2)!y₀^{1−k}] and Borel-resum it exactly:
Σ_{k≥K}2Re[(k−2)!y₀^{1−k}]h^{k+1} → R^{lead}_K = 2Re[h² i^{1−K}∫₀^∞e^{−πpσ}σ^{K−2}(1+iσ)^{−1}dσ]. Laplace's method at
σ* = (K−2)/(πp) ≈ 1 gives R^{lead}_K ≈ 2Re[i^{1−K}/(1+i)]h²Γ(K−1)(πp)^{1−K} = √2 cos(π(K−1)/2 + π/4)·Γ(K−1)/(p²(πp)^{K−1}),
and |√2 cos(π(K−1)/2 + π/4)| = 1 for every integer K. With K − 1 = πp + ν, Stirling gives
  |R_K(p)| ≈ √2 p^{−5/2} e^{−πp} e^{ν²/(2πp)}.
The second large-order term (the log k term, relative size O(log k/k)) changes this by a factor 1 + O(log p/p) whose sign
depends on the parity of K; it is not small at moderate p. *Numerical* (Table C2, p = 4, 4.5, …, 14): |R_{K_s}|/ref =
0.62–0.86 with K_s always even, slowly increasing with p (0.70 at p = 6.5, 0.77 at 9, 0.79 at 11, 0.82 at 13.5), consistent
with 1 − O(log p/p); at p = 3 and 3.5, K_s = 8 is far from πp + 1 and the ratio is 11 and 17.5. At odd K (Table C) R_K/R^{lead}_K = 1.41–1.67, at even K 0.60–0.85; both approach 1 slowly. Unproved step:
the uniform control of the remainder of the large-order expansion inside the Borel integral (Dingle–Berry type terminant
analysis), which B2 does with a non-sharp constant only.

So the "true constant and power": |error at the smallest term| = √2·p^{−5/2}·e^{−πp}·(1 + O(log p/p)) (DERIVED), about half
the smallest odd term 2√2 p^{−5/2}e^{−πp}. The rigorous bound B2 has the same exponential rate e^{−πp} and the power
p^{−5/2} up to a factor O(log p).

---

## 5. Task C — the first beyond-all-orders term — DERIVED + NUMERICAL

There is no Stokes line on the positive axis (B1), so the "exponentially small correction" is the optimally truncated remainder
itself. Its leading term:

  R_K(p) = √2 cos(π(K−1)/2 + π/4) · Γ(K−1) / (p² (πp)^{K−1}) · (1 + O(log p/p)),  K − 1 = πp + O(√p).

The O(log p/p) correction depends on the parity of K and is still 15–50 % at p ≈ 10 (Table C); the refined model below
captures it. The term oscillates with K with period four (signs + − − + for K ≡ 1, 2, 3, 0 mod 4) and has the Gaussian envelope
√2p^{−5/2}e^{−πp}e^{(K−1−πp)²/(2πp)}. As p grows with K = K_s(p) (always even), the sign follows K mod 4, e.g. K_s = 28, 30,
32, 34, 36 at p = 8, 9, 9.5, 10, 10.5 gives signs +, −, +, −, + (Table C2: all 23 values of p from 3 to 14 have the predicted
sign).

**Refined model (C2).** Replace the tail model by the exact leading singular part of the n = 1 symbol,
A ≈ (y₀ − y)^{−1−h} + c.c. (all orders in h): M(h) = 2Re[h²∫₀^∞e^{−v/h}(iπ − v)^{−1−h}dv] (an incomplete Gamma function),
whose expansion coefficients S_k are the l = 0 singular coefficients (S_k/N_k = 0.939, 0.988, 0.995 at k = 20, 40, 60). Then
R^{sing}_K := M(h) − Σ_{k<K}S_kh^{k+1} contains the log-term (parity) correction automatically. *Numerical*: R_K/R^{sing}_K =
1.017–1.040 (p = 8), 1.012–1.023 (p = 10), 1.009–1.015 (p = 12) for nine consecutive K around πp at each p; the deviation
is O(1/p), as the next singular term (l = 1, relative size 1/k) predicts. DERIVED (the l ≥ 1 and n ≥ 2 contributions are
bounded heuristically, not proved).

**Table C** (p = 10; R_K = actual remainder):

| K | R_K | R/R^{lead} | R/R^{simple} | R/R^{sing} |
|---|---|---|---|---|
| 30 | −8.443e−17 | 0.740 | 0.724 | 1.0159 |
| 31 | −1.720e−16 | 1.506 | 1.598 | 1.0200 |
| 32 | 7.864e−17 | 0.777 | 0.765 | 1.0145 |
| 33 | 1.514e−16 | 1.497 | 1.492 | 1.0175 |
| 34 | −8.226e−17 | 0.808 | 0.796 | 1.0132 |
| 35 | −1.515e−16 | 1.489 | 1.396 | 1.0154 |

(R^{simple}_K = √2cos(π(K−1)/2+π/4)Γ(K−1)/(p²(πp)^{K−1}).) Rows for p = 8, 12 and Table C2 are in ru3_log.txt.

---

## 6. What the paper may say (statement-ready)

For other rulers the answer depends only on whether $q$ stays bounded. If $q=q(1/p)$ is bounded and analytic in $1/p$, the
series converges; for the dual ruler it is even exact: $L_{p/(p-1)}(p)=\frac4pB(\frac1p,\frac1p)=8\,\Gamma(1+\frac1p)^2/\Gamma(1+\frac2p)$,
twice the area of the $p$-disc, so $L_{p/(p-1)}(p)/8=1-\zeta(2)p^{-2}+2\zeta(3)p^{-3}-\dots$ converges exactly for $p>1$. If $q$ grows
linearly, $q=\lambda p+c+O(1/p)$ with $\lambda>0$ and analytic in $1/p$, the series diverges for every $p$:
$N_k=-\frac{2(k-2)!}{\pi^{k-1}\lambda^k}\bigl[\cos\bigl(\pi(\lambda-c)-\tfrac{\pi(k-1)}2\bigr)+O\bigl(\tfrac{\log k}k\bigr)\bigr]$, and the Borel
singularities sit at $\pm i\pi\lambda$, that is, at $\pm i\pi$ in the variable conjugate to the ruler $q$. For $q=p$ the series is
Borel summable, its Borel sum is $\pi_p/4$ for $p>(1+\sqrt5)/2$, and truncating at the smallest term leaves an error of size about
$\sqrt2\,p^{-5/2}e^{-\pi p}$ (we prove a bound of order $e^{-\pi p}$ times a power of $p$; the constant $\sqrt2$ is derived and checked numerically).

---

## 7. Remarks, limits, and what is not done

- A3 is proved by adapting PIPDIV (refereed); the O-constants are not computed, as there.
- In A3 (ii), for non-linear Q, the real form needs the complex b₁ term; the complex form displayed first is the statement.
- B1's threshold p > φ comes from the crude bound in Lemma P (Cauchy radius ρ < 1); we expect the Borel sum to equal π_p/4 − 1
  for every p > 1 (the true exponential type of 𝔅 on [0,∞) should be 1), not proved.
- B2 is rigorous but not sharp (factor ≈ log p and constants); the sharp constant √2 is DERIVED (B3).
- A4 (d ≥ 2) and the non-analytic rulers (√p, log p) are DERIVED only.

## 8. Scripts and logs (run with ../../../.venv/bin/python, from this folder)

- ru0_Jtable.py K DPS DEG S — J_{i,j} table; ru0_J_K90_dps100_deg6_S700.txt (log ru0_log_a.txt) and
  ru0_J_K90_dps80_deg5_S640.txt (ru0_log_b.txt); J_{1,j} checked against j!(1−2^{−j−1})ζ(j+2) (9.6·10⁻¹⁰⁰, 9.3·10⁻⁸⁰).
- ru1_family.py FAM JFILE K DPS — N_k for q = λp (lam:λ), λp + c (lin:λ:c), p² (sq); exact n = 1 term T_k, singular part S_k,
  two-term law; files ru1_N_*_K90.txt (run a with the 100-digit table for λ = 1, 2, ½ and q = p+1; run b for λ = 3/2,
  q = p+½, q = p²), runb/ (all run-b files); logs ru1_log_*_{a,b}.txt; ru1_crosscheck_log.txt.
- ru3b_p20_log.txt — p = 20: π_20/4, K_s = 66, R_{K_s} = −3.637·10⁻³¹ (ref 4.078·10⁻³¹, ratio 0.89), R/R^{sing} = 1.0035.
- ru2_direct.py — quadrature checks, dual ruler, Table A′; ru2_log.txt.
- ru3_borel.py — π_p values, truncation tables, models R^{lead}, R^{simple}, R^{sing}, Tables B, C, C2 (its bound column is the
  crude Lemma-R corollary); ru3_log.txt.
- ru4_bound.py [p …] — Theorem B2 evaluated; ru4_log.txt (p = 5, 10, 20), ru4_log_p3.txt, ru4_log_p2.txt (p = 2, small
  admissible grid run inline with the functions of ru4_bound.py).
- ru5_tables.py — Table A; ru5_log.txt.

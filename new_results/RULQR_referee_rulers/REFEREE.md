# RULQR — referee report on RULQ (rulers that grow with p; Borel summability of π_p)

Refereed: `../RULQ_growing_rulers/result.md` and its scripts. Background read: `paper/ieee/main.tex` (thm:pip, prop:ruler
and the text around them), `../PIPDIV_divergence_proof/result.md`, `../PIPDIVR_referee_divergence/REFEREE.md`.
Notation as in RULQ: h = 1/p, L_q(p)/8 ~ 1 + Σ_k N_k h^{k+1}, y₀ = iπ. All my numbers come from the scripts in this folder
(Sec. 7). I did not use RULQ's J-tables or series code.

## Verdicts

| Claim | RULQ label | Verdict |
|---|---|---|
| A2 bounded analytic rulers converge (Morera) | PROVED | **CONFIRMED** |
| A1 dual ruler L_{p′}(p) = (4/p)B(1/p,1/p) = 8Γ(1+h)²/Γ(1+2h) = 2·Area; series, radius 1 | PROVED | **CONFIRMED** |
| A3 linear growth q = pQ(1/p): divergence for every p, two-term law, Borel radius πλ, singularities ±iπλ | PROVED | **CONFIRMED** (three wording repairs, R3a–c) |
| "converges iff q stays bounded" (RULQ §0, §6 first sentence) | PROVED for d = 0, 1 | **NOT PROVED** as stated in §6; proved only for rulers q = p^d Q(1/p), d ∈ {0, 1} (repair R4) |
| A4 q = p² diverges, Gevrey order 1/2 | DERIVED | **NOT PROVED** (label DERIVED is right; numerics reproduced) |
| B1 Borel transform: radius π, strip |Im t| < π, Laplace sum = π_p/4 − 1 for p > φ | PROVED | **CONFIRMED** (threshold φ correct as proved; not sharp) |
| B2 explicit truncation bound | PROVED | **CONFIRMED AFTER REPAIR** (R5: range of δ₂; status of the numbers) |
| B3 true error ≈ √2 p^{−5/2}e^{−πp}; C sign pattern of period 4 | DERIVED | **NOT PROVED** (labels right; numerics reproduced) |

## 1. Independent numerical method

From the definition L_q(p) = (4/p)∫₀¹(t^{q/p−q} + (1−t)^{q/p−q})^{1/q}dt, symmetry, t^{−a} factored out on [0,½]
(a = r = q(1−h)), u = t/(1−t) and u = e^{−w/r}, I get, for p > 1, q > 0, η = 1/q:

  L_q(p)/8 = 2^{−h} + Ê(h),  Ê(h) = (h/r)∫₀^∞ e^{−hw/r}(1+e^{−w/r})^{−1−h}[(1+e^{−w})^η − 1]dw.

(This is RULQ's (2.1)/(2.3); the derivation is forced.) What is independent is the extraction: I evaluate G(h) = Ê(h)/h³
(Ê/h for a bounded ruler, Ê/h⁵ for q = p²) in **Arb** (`acb.integral`, certified error radius ≈ 10⁻²⁴⁵) at Chebyshev nodes on
[0, H], and read the Taylor coefficients off the Chebyshev interpolant. Since the Borel singularities are at ±iπλ, the
non-perturbative part is O(e^{−πλ/H}) = O(e^{−785}) for H = 0.002, λ = ½, far below the amplified rounding. Two runs
(H = 0.002, 110 nodes, 220 digits; H = 0.0015, 130 nodes, 260 digits) agree to all 60 printed digits for k ≤ 45.

Results (`rr3_log.txt`, `rr3b_log.txt`, `rr3c_log.txt`):
- λ = 2, λ = ½, q = p + 1, q = p, q = p²: N_k agree with RULQ's `ru1_N_*` files to ≤ 2·10⁻⁵⁹ relative (k ≤ 45; for q = p², k ≤ 41).
  For q = p they agree with PIPDIVR's ball-arithmetic N_k to all printed digits.
- Two-term law of A3(ii): r_k = (ν_k − F_k)k²/log²k stays bounded for k ≤ 50 in every family.
- **Generic phase.** All of RULQ's test families (λ = 2, ½, 3/2, q = p+1, p+½, p) have 2(λ − c) ∈ ℤ, where one parity is
  suppressed. I added q = p + 1/3 + 1/p (λ = 1, c = ⅓, b₁ = Q″(0)/2 − c = ⅔), where both parities lead with weights
  cos = 0.866 and 0.5, and the b₁ term is non-zero: ν_k → F_k (k = 48: 0.990 vs 0.983; k = 49: 0.841 vs 0.793), r_k bounded
  (≈ 1–8). This tests the phase e^{iπ(1+λ−c)} and the b₁y₀ term, which RULQ's data did not.
- Bounded ruler q = 2 + 1/p: |N_k|^{1/k} ≈ 0.67–0.76 for k ≤ 50 (geometric, as A2 says).
- q = p²: |N_k|^{2/k}/(k/(2eπ)) = 1.06–1.27 for 20 ≤ k ≤ 60, as RULQ reports.

Direct quadrature (`rr4_log.txt` (a)): the geometric definition with x = cos^{2/p}t, y = sin^{2/p}t, split at π/4
(t = u^{p/2} removes the t^{2/p−1} endpoint), and a graph form x = (1−y^p)^{1/p}, agree with the Arb w-form to 40 digits at
(p, q) = (3,2), (5,5), (4,8), (10,5), (6,7), (2.5,2.4), (20,40), and with 8Γ(1+1/p)²/Γ(1+2/p) at (4,4/3), (1.5,3), (7,7/6).
Table A′-type check (b): |L/8 − sum up to the smallest term| is of the size of the smallest term in all 8 cases
(λ = 2 at p = 4, 6; λ = ½ at p = 8, 16; q = p+1 at p = 8, 12; generic at p = 8, 12).

## 2. A1 (dual ruler) and A2

**A1.** RULQ's proof is correct: q = p′ gives q/p − q = −1 and 1/q = 1 − h, so the integrand is (t(1−t))^{h−1} and
L = (4/p)B(h,h) = 8Γ(1+h)²/Γ(1+2h). The series: log Γ(1+h)²/Γ(1+2h) = Σ_{k≥2}(−1)^k(2 − 2^k)ζ(k)h^k/k gives
1 − ζ(2)h² + 2ζ(3)h³ + (−7π⁴/180 + π⁴/72)h⁴ = … − (π⁴/40)h⁴ (hand check; `rr6_log.txt` checks D₂, D₃, D₄ to 40 digits).
Radius: Γ(1+h)² has double poles at h = −1, −2, …; 1/Γ(1+2h) has simple zeros at h = −½, −1, −3/2, …; so h = −½ is a zero
(not a singularity) and h = −1 is a simple pole, residue −2 (`rr6`: (−1)^k D_k → −2). Radius exactly 1; the series converges
exactly for p > 1. CONFIRMED.

**The planner's one-line proof is right.** Let K be the p-disc, s Euclidean arclength, T the unit tangent, n the outward unit
normal (n is T turned by 90°). Then ∫_{∂K} h_K(n) ds = ∫_{∂K} x·n ds = ∫_K div x dA = 2·Area(K); h_K(n) = ‖n‖_{p′}; and
‖·‖_{p′} is invariant under (a, b) ↦ (−b, a), so ‖n‖_{p′} = ‖T‖_{p′} and the left side is L_{p′}(p). (It holds for every
1 ≤ p ≤ ∞; at p = ∞ it is L₁(∞) = 8, consistent with prop:ruler.) RULQ's Beta-integral proof is a different route to the same
identity; the two agree (`rr6`: ∫h_K(n)ds = 2·Area to 25 digits at p = 3). In Minkowski geometry this is the standard
mixed-area formula for Minkowski length (length of ∂C in the norm with unit ball B equals 2V(C, B° turned by 90°)):
**classical; no reference verified.**

**A2.** Checked line by line. For h in U: |u^{r(h)}| = u^{Re r} < 1 on (0,1), so for h in a compact K ⊂ U, z stays in a compact subset
of the closed unit disc that avoids −1; |Log(1+z)| ≤ C_K|z|; |Ψ| ≤ C′_K u^{Re r}; integrand ≤ C″_K u^{Re h + Re r − 1}, uniformly on compact
K ⊂ U; Morera. The fixed-q and dual cases of U are right ({Re h < 1} resp. {|h| < 1}). Only Q(0) > 0 is used (Q(0) ≥ 1 is
there so that q is a norm). CONFIRMED.

## 3. A3 (linear growth), step by step

I re-derived each adaptation of PIPDIV.
1. Laplace form: r = λ/h + β(h), β(0) = c − λ = b₀, β′(0) = Q″(0)/2 − c = b₁; e^{−nrv} = e^{−nλv/h}e^{−nβv}, so
   H_n = e^{(1−nβ)v}(1+e^v)^{−1−h}; for Q ≡ 1 this is PIPDIV's H_n. Σ|C(η,n)| = 1 for 0 < η ≤ 1. Correct.
2. **n = 1 dominance.** With t = nλv the symbol is e^{t/(nλ)}e^{−β(h)t/λ}(1+e^{t/(nλ)})^{−1−h}; the factor e^{−βt/λ} does not
   depend on n, so the bounds are uniform in n. Singularities at t = ±iπλn(2j+1): radius πλ for n = 1, 2πλ for n ≥ 2. The tail
   bound e^{Bt}, B = sup|β|/λ, is right. U_k = O(C^k + k(k−2)!(1.9πλ)^{−k}) — *cosmetic repair:* "2^k" in (i) should be C^k
   (the h-disc where Q, β, η are controlled has radius h₁, not ½). Relative to (k−2)!(πλ)^{−k} this is still geometric.
3. **Branch points.** 1 + e^v = (y₀ − v)φ(v − y₀) because e^{y₀} = −1; the extra factor e^{(1−β)v} is entire, so the branch
   bookkeeping of PIPDIV (with the PIPDIVR repair R1, domain W) is unchanged. Singular points ±iπ in v, i.e. ±iπλ in t.
4. **Phase.** Residue of e^{(1−b₀)v}/(1+e^v) at iπ is −e^{(1−b₀)y₀}, matching the l = 0 term; e^{(1−b₀)y₀} = e^{iπ(1+λ−c)}.
   For q = p this is e^{2πi} = 1 (PIPDIV). Checked numerically in five families including the generic one.
5. **Two-term law.** m![v^m] of e^{(1−β)y₀}(y₀−v)^{−1−h} = e^{(1−β)y₀}y₀^{−1−m}e^{−hLog y₀}(1+h)_m; with the weights λ^{−m−1}
   and 1/Q(h): [h⁰] gives e^{(1−b₀)y₀}λ^{−k}(k−2)!y₀^{1−k}; [h¹] adds (k−3)!λ^{1−k}y₀^{2−k}(ψ(k−2)+γ−Log y₀ − b₁y₀ − c/λ).
   Real form and the λ = 1, c = 0 check (1 + π²/(2(k−2)) at odd k) are right.
6. **Uniformity.** Cauchy on |h| = 1/log k needs only analyticity of Q, 1/Q, β on a fixed disc; fine. But the hidden constants
   grow with λ: r_k ≈ 20 for λ = 2, ≈ 5 for λ = 1, ≈ 1 for λ = ½ (roughly ∝ λ², the next term being O((πλ log k/k)²)).
   The law is accurate only for k ≫ πλ log k (R3c).
7. **(iii) divergence for every p.** cos(θ − π(k−1)/2) and cos(θ − πk/2) differ by a quarter turn, so one of each consecutive
   pair is ≥ 1/√2 in modulus. **This also covers the degenerate case.** When 2(λ − c) ∈ ℤ (e.g. λ = 1, c = 0; every family in
   RULQ's Table A), the leading cosine is 0 for one parity, and for that parity N_k is the second term, of order
   (k−3)!log k/(π^{k−2}λ^{k−1}); the other parity has |cos| = 1 and gives divergence. So "diverges for every p" is proved in
   all cases. The §6 formula "−2(k−2)!/(π^{k−1}λ^k)[cos(…) + O(log k/k)]" is true in the degenerate case too, but for the
   suppressed parity it says only "O((k−3)!log k …)" — it should say so explicitly (R3a).
8. **(iv) Borel radius and singularities.** |N_k|/k! ≍ k^{−2}(πλ)^{−k} along a subsequence, ≲ k^{−2}(πλ)^{−k}log k always:
   radius exactly πλ (normalisation Σ N_k t^k/k! ↔ Σ N_k h^{k+1} via ∫e^{−t/h}t^k dt = k!h^{k+1}). The location (only ±iπλ on
   the circle) uses the strip analyticity of the remark A3′, which RULQ only sketches; I checked that the sketch works:
   |H_n(t/(nλ),η)| ≤ e^{πρ}σ^{−1−ρ}e^{(B + ρ/λ)|t|} on |Im t| ≤ λ(π − δ), uniformly in n, and Lemma P goes through. The
   constant in "p larger than an explicit constant" is not computed; "explicit" should go. "The Borel singularities sit at
   ±iπλ" should read "the nearest singularities" (there are others at ±iπλ(2j+1) and from n ≥ 2) (R3b).

**R4 (scope).** RULQ §3.4 states the scope correctly (proved for d = 0, 1; derived for d ≥ 2; non-analytic rulers outside).
The paper text in §6 ("the answer depends only on whether q stays bounded") drops all of that. It is not proved for
q = p² (A4 is derived), not for rulers that are not of the form p^dQ(1/p) (q = √p, log p, or bounded rulers like 2 + p^{−1/2},
which have no 1/p power series at all), and the paper's §6 sentence should be restricted to analytic rulers.

## 4. B1 (Borel summability)

Checked every inequality of Steps 1–5: |1+e^v|² ≥ (1 − cos δ)(1+e^{2x}) (difference cos δ(1−e^x)² ≥ 0 for δ ≤ π/2);
|(1+e^v)^{−1−η}| ≤ σ^{−1−ρ}e^{πρ}(1+e^{2x})^{−(1−ρ)/2}; |H_n(t/n,η)| ≤ e^{πρ}σ^{−1−ρ}e^{(1+ρ)(Re t)⁺}; |C(η,n)| ≤ ρe^ρ n^{ρ−1};
Lemma P (|B̂| ≤ C(1+|t|/ρ)e^{max(1+ρ,1/ρ)|t|}); identification by Watson + uniqueness; β₀ part entire. Normalisation is right.
**Threshold.** The Laplace integral needs p > max(1+ρ, 1/ρ) for some ρ ∈ (0,1); min_ρ max(1+ρ,1/ρ) = φ at ρ = 1/φ. So
"every p > (1+√5)/2" is exactly what the proof gives. It is not sharp: for the n = 1 part the Borel transform is
B₁(t) = ∫₀^t g(y)√((t−y)/ℓ(y)) J₁(2√((t−y)ℓ(y))) dy (g = e^{2y}/(1+e^y), ℓ = log(1+e^y) ≥ log 2, |J₁| ≤ 1), of exponential type 1
on [0,∞); if the whole 𝔅 has type 1 + ε there, analytic continuation in p gives every p > 1 (RULQ's expectation). Not proved.
CONFIRMED.

## 5. B2 (explicit bound)

**Derivation.** Lemma R checked (Cauchy remainder on the rectangle; left edge → 0; horizontal edges |s − y| ≥ r and
|f(x − ir)| = |f(x + ir)| by reality; right edge; Tonelli). Theorem B2: the n = 1 part (g, ℓ bounds, Cauchy in η, m₁ exact),
the index shifts N = K−2−j and N = K−j, the N ≤ 1 geometric tails, R0 — all correct.
**Repair R5.** The n ≥ 2 part applies Step 2(a) of B1 with δ = δ₂/2, which needs δ₂/2 ≤ π/2. The stated range
0 < δ₂ < 2π is wrong; it must be 0 < δ₂ ≤ π (for δ ∈ (π/2, π) the inequality |1+e^v|² ≥ (1 − cos δ)(1+e^{2x}) fails). All
evaluated choices (δ₂ ∈ {1, 1.5, 2}) satisfy it. Typo "0 < r < π′" in Lemma R.

**My evaluation** (`rr5_bound.py`, own code): same parameters, two ways on the same grids — midpoint sums (ordinary numerics)
and **rigorous upper Riemann sums in Arb**: on each cell, Γ(N+1,·) and (x²+r²)^{−N/2} are evaluated at the cell end that
maximises them (both monotone), m₁ is evaluated on the cell as a ball, all products rounded up; closed-form tails for
x < −50r, x > X_R and y > Y (constants re-derived; I keep the 1/(πr) factor RULQ's code drops, which only made theirs larger).

| p | K | (ρ, δ, ρ₂, δ₂) | RULQ (numerical) | mine, midpoint | **mine, certified upper sum** | true |R_K| | ratio |
|---|---|---|---|---|---|---|---|
| 5 | 18 | 0.6, 0.5, 0.6, 2 | 1.518·10⁻⁶ | 1.518·10⁻⁶ | **1.583·10⁻⁶** | 2.626·10⁻⁹ | 603 |
| 10 | 34 | 0.35, 0.05, 0.6, 1 | 1.596·10⁻¹⁵ | 1.596·10⁻¹⁵ | **1.613·10⁻¹⁵** | 8.226·10⁻¹⁷ | 19.6 |
| 20 | 66 | 0.2, 0.05, 0.9, 1.5 | 5.407·10⁻³⁰ | 5.407·10⁻³⁰ | **5.454·10⁻³⁰** | 3.637·10⁻³¹ | 15.0 |

(RULQ did not print ρ₂, δ₂ for p = 5, 10; I searched the same grid and reproduce their totals and R2 parts.)
**Is "proved" fair?** The inequality R0 + R1 + R2 is proved. RULQ's numbers were not certified, so in RULQ the numerical values
are "a proved bound, evaluated numerically". With the upper sums above the values are certified (modulo correctness of my
~150-line script, which reproduces RULQ's numbers in midpoint mode). Also: that the bound is O(p^{−5/2}log p·e^{−πp}) for all
large p is DERIVED in RULQ, not proved; the paper must not say "we prove a bound of order e^{−πp} times a power of p" without
either that asymptotic analysis or the qualifier "at p = 5, 10, 20". Suggested wording: "We prove an explicit bound for the
truncation error. Evaluated in ball arithmetic at the smallest term, it gives 1.6·10⁻⁶, 1.6·10⁻¹⁵ and 5.5·10⁻³⁰ at
p = 5, 10, 20, within factors 600, 20 and 15 of the true error."

## 6. B3 / C (derived) — reproduced

π_p/4 from Arb (60 digits) at p = 5, 10, 20 agrees with RULQ to all 40 printed digits. Optimal truncation (`rr4_log.txt` (c)):
p = 5: K_s = 18, R = −2.626·10⁻⁹, K_opt = 14, R = −1.620·10⁻⁹; p = 10: K_s = 34, R = −8.226·10⁻¹⁷, K_opt = 32,
R = 7.864·10⁻¹⁷; p = 20 (PIPDIVR N_k): K_s = 66, R = −3.637·10⁻³¹. |R_{K_s}|/(√2p^{−5/2}e^{−πp}) = 0.69, 0.81, 0.89.
R_K/R^{simple}_K is positive for all 27 K checked (so the predicted signs + − − + for K ≡ 1, 2, 3, 0 mod 4 hold), and
alternates between ≈ 0.3–0.9 and ≈ 1.2–2.6 with parity, approaching 1 slowly (0.87/1.21 at p = 20). DERIVED labels are right.

## 7. Paper text (RULQ §6) — sentence by sentence

1. "For other rulers the answer depends only on whether q stays bounded." — **Overclaim** (R4). Proved only for rulers analytic
   in 1/p with d ∈ {0,1}.
2. "If q = q(1/p) is bounded and analytic in 1/p, the series converges" — proved (A2). Add q ≥ 1 (or Q(0) > 0).
3. Dual ruler formula, "twice the area", series, "converges exactly for p > 1" — proved (A1).
4. "If q grows linearly … the series diverges for every p" — proved (A3).
5. N_k formula — proved; must add that it is for large k and that for 2(λ − c) ∈ ℤ one parity is smaller by log k/k (R3a).
6. "the Borel singularities sit at ±iπλ" — say "nearest"; the "variable conjugate to q" gloss is interpretation, keep it out
   or mark it as such.
7. "For q = p the series is Borel summable, its Borel sum is π_p/4 for p > (1+√5)/2" — proved (the threshold is the proof's).
8. "truncating … error of size about √2 p^{−5/2}e^{−πp}" — derived, not proved; "(we prove a bound of order e^{−πp} times a power
   of p …)" — overclaim (Sec. 5).

**Rewritten text (statement-ready):**

> Which rulers give a convergent series? We can answer this for rulers that depend analytically on $1/p$. Write
> $L_q(p)/8\sim1+\sum_k N_kp^{-k-1}$. If $q=Q(1/p)$ with $Q$ analytic at $0$ and $Q(0)\ge1$, the series converges; the proof
> of Proposition~\ref{prop:ruler} still works when $q$ depends analytically on $1/p$. For the dual ruler $q=p/(p-1)$ the length
> is even exact: $L_{p/(p-1)}(p)=\frac4pB(\frac1p,\frac1p)=8\,\Gamma(1+\frac1p)^2/\Gamma(1+\frac2p)$, twice the area of the
> $p$-disc. (The $\ell^{p/(p-1)}$ norm is the support function of the $p$-disc and does not change when we turn a vector by
> $90^\circ$, so this length is $\int h_K(n)\,ds=2\,\mathrm{Area}(K)$.) Hence
> $L_{p/(p-1)}(p)/8=1-\zeta(2)p^{-2}+2\zeta(3)p^{-3}-\frac{\pi^4}{40}p^{-4}+\dots$, and this series converges exactly for $p>1$.
> If the ruler grows linearly, $q=\lambda p+c+O(1/p)$ with $\lambda>0$ and analytic in $1/p$, the series diverges for every
> $p$. For large $k$,
> $N_k=-\frac{2(k-2)!}{\pi^{k-1}\lambda^k}\bigl[\cos\bigl(\pi(\lambda-c)-\tfrac{\pi(k-1)}2\bigr)+O\bigl(\tfrac{\log k}k\bigr)\bigr]$.
> When $2(\lambda-c)$ is an integer, as for $q=p$, the cosine vanishes for one parity of $k$, and those coefficients are
> smaller by a factor of order $\log k/k$. The Borel transform $\sum_kN_kt^k/k!$ has radius of convergence $\pi\lambda$, and
> its only singularities on that circle are $\pm i\pi\lambda$. For faster growth, such as $q=p^2$, we expect the series to
> diverge too, with $N_k$ growing like $(k/2)!$; we have derived this and checked it numerically, but not proved it.
> For $q=p$ the series is Borel summable: its Borel sum is $\pi_p/4$ for every $p>(1+\sqrt5)/2$ (we expect every $p>1$).
> We also prove an explicit bound for the error of the truncated series. At the smallest term, evaluated in ball arithmetic,
> it gives $1.6\cdot10^{-6}$, $1.6\cdot10^{-15}$ and $5.5\cdot10^{-30}$ at $p=5,10,20$, within factors $600$, $20$ and $15$ of
> the true error. The true error there is about $\sqrt2\,p^{-5/2}e^{-\pi p}$; this is derived, not proved, and we checked it
> numerically for $4\le p\le20$.

(If the paper does not want to rely on RULQR's ball-arithmetic evaluation, replace "evaluated in ball arithmetic" by
"evaluated numerically".) The last sentence of the current paper paragraph ("We do not know what happens for other rulers that
grow with p") must then go, and the summary table needs rows for A1/A2 (proved), A3 (proved), B1 (proved), B2 (proved;
values certified by RULQR), A4/B3 (derived).

## 8. Scripts and logs (run from this folder with `../../../.venv/bin/python`)

- `rr1_values.py FAM H D DPS OUT` — Arb values of G = Ê/h^s at Chebyshev nodes (families lam2, lamhalf, p1, p, gen, bdd, sq);
  data `rr1_G_*_{a,b}.txt` (a: H = 0.002, 110 nodes, 220 digits; b: H = 0.0015, 130 nodes, 260 digits).
- `rr2_extract.py VALFILE H KMAX OUT [s]` — Chebyshev → Taylor → N_k; `rr2_N_*_{a,b}.txt`.
- `rr3_compare.py` → `rr3_log.txt` (a vs b, vs RULQ, vs PIPDIVR, two-term law); `rr3b_generic.py` → `rr3b_log.txt`
  (generic phase, bounded ruler); `rr3c_square.py` → `rr3c_log.txt` (q = p²).
- `rr4_direct.py` → `rr4_log.txt` (direct quadratures, Table A′-type check, π_p and optimal truncation at p = 5, 10, 20).
- `rr5_bound.py p K ρ δ ρ₂ δ₂ c` → `rr5_log_p{5,10,20}.txt` (Theorem B2, midpoint and certified upper sums).
- `rr6_dual.py` → `rr6_log.txt` (dual-ruler coefficients, residue, support-function identity).

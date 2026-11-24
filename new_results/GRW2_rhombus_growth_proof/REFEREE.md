# Referee report on GRW2 (`proof.md`): radius zero at the diamond end, every rhombus

Target: `paper/ieee/main.tex`, Theorem "Radius zero" (`thm:zero`, line 1229) and its proof sketch.
Background taken as correct: `GRW1_diamond_growth_borel/proof.md` (A = B).
The referee's own scripts are in `referee/` (run `ref_coeffs.py` first, which writes `In_direct.pkl`, then `ref_analysis.py` and
`ref_square_next_term.py`; output in `referee/ref_output.txt`). The author's `numerics.py` was not run; its output
file was read only after the referee's numbers were in.

Summary: no step of the proof of Theorem A (the paper's theorem) fails. Two places need small wording fixes (items 6a, 6b).

---

## 1. Integral representation and Borel set-up: HOLDS

Re-derived from the definition. The arc x = A t^{1/p}, y = B(1−t)^{1/p} satisfies |x/A|^p+|y/B|^p = 1 and has
speed (1+η)√(A²t^{2η}+B²(1−t)^{2η}), so I(η) = (1/r)∫₀¹√(A²t^{2η}+B²(1−t)^{2η}) dt. With t/(1−t) = e^l,
t^{2η} = w^η e^{ηl}, (1−t)^{2η} = w^η e^{−ηl}, dt = w dl, and folding l ↦ −l gives exactly
I = ∫₀^∞ e^{−u}H(ηu)K(η,u) du with h₁(z) = √(A²+B²e^{−2z})/r, h₂(z) = √(B²+A²e^{−2z})/r. The identification with
the paper's H = e^{−z/2}[g(z)+g(−z)] is correct: (1±q)/2 = A²/r², B²/r², so r²(cosh z+q sinh z) = A²e^z+B²e^{−z}.
The integrability bound (H(x) ≤ (A+B)(1+e^{−x})/r, K ≤ 1 for η > −2) gives −1 < η < 1.

The coefficient formula (Lemma 2) is the Leibniz rule at η = 0: ∂ⁿ_ηF(0,u) = n!Σ_j H_{n−j}u^{n−j}k_j(u). The
dominated differentiation is correctly justified (Cauchy estimates for H on disks of radius π/4 in the strip,
|∂^j_ηK| ≤ log^j 2 ≤ 1). M_{n,0} = n!η_D(n) and the bound |M_{n,j}| ≤ (n−j)!/(j!(j+1)^{n−j+1}) are right.
The Borel transform is used only as the formal series Σ(I_n/n!)ζⁿ, which is all Cauchy–Hadamard needs.

Numerical: the t-form and the Laplace form agree to 20 printed digits at η = 0.3, −0.3, 0.7 for all three
rhombi. I₀ = 1 exactly.

## 2. Singularity analysis: HOLDS

- Zeros of R₁: e^{−2z} = −A²/B² gives z = L + iπ(k+½), all simple (R₁′ = −2B²e^{−2z} ≠ 0); R₂: −L + iπ(k+½).
  The smallest modulus is ρ, reached only at k ∈ {0,−1}, and the next is ρ₂ = √(L²+9π²/4). Checked.
- No singularity closer than ρ (the upper bound): on D(0,ρ), which is simply connected, R₁ and R₂ have no zeros,
  so h₁ and h₂ are analytic there. They agree with the principal roots on D(0,ρ₂) ∩ S, which is convex because the
  slits t·z_j (t ≥ 1) have |Im| ≥ π/2. The slit disk Ω is star-shaped as claimed. Correct.
- Each z_j is a genuine branch point (Lemma 3(e)): if H extended to a disk at z₁, then h₁ = H − h₂ would extend
  (h₂ is analytic near z₁ when A ≠ B; D(0,ρ) ∩ D(z₁,ε) is convex, so the continuations match), h₁² = R₁/r²
  would persist, and a simple zero cannot be a square. A = B: H = 2h₁. Correct.
- Local constants: c_j as stated. The referee checked them by evaluating h_i(t z_j)/√(1−t) at t = 1 − 10⁻¹⁴ for all
  four points and all three rhombi. All agree to 12 digits, including the branch choice. The expansion
  ψ₁ = 1 + (z₁/2)s + O(s²) follows from φ(w) = (1−e^{−2w})/(2w) = 1 − w + O(w²).
- The "K-part" (Lemma 4): the exact claim is |E_n| ≤ C′_σ(2σ)^{−n} for all n ≥ 1, every σ < ρ, with
  C′_σ = C_σ[1 + ½Σ_{j≥1}(2σ)^j/(j!)²] and C_σ = max_{|z|=σ}|H|. The constant is uniform in n. Every inequality checks:
  |η_D(n)−1| ≤ 2^{−n} (alternating, decreasing); n!/(n−j)! ≥ j!; (j+1)^{−(n−j+1)} ≤ 2^{−(n−j+1)} for j ≥ 1. So
  𝓔 has radius ≥ 2ρ. The referee also checked the inequality with explicit numbers (σ = 0.9ρ, C_σ by sampling 720
  points): it holds for every n ≤ 40 in all three cases. E_n(2ρ)ⁿ tends to 0 slowly, like n^{−3/2}, as
  E_n ≈ −2^{−n}H_n predicts, so 2ρ is the true radius of 𝓔.
- Hidden interchanges: there are none. The only termwise integration is M_{n,0} (Tonelli for n ≥ 2, direct for
  n = 0, 1). The K-expansion is used only through derivatives at η = 0. Lemma 4 bounds a finite sum by a convergent
  series. 𝓑 = H + 𝓔 is an identity of formal power series.

## 3. From singularities to limsup = 1/ρ, and the cancellation question: HOLDS

Cauchy–Hadamard is applied to the right series, 𝓑(ζ) = Σ(I_n/n!)ζⁿ = H + 𝓔 coefficientwise. R(H) = ρ exactly
and R(𝓔) ≥ 2ρ, so R(𝓑) = ρ, and limsup|I_n/n!|^{1/n} = 1/ρ is an identity, not an estimate. Interference among the
four points (or the two merged points when A = B) cannot affect it. The lower bound comes entirely from the
non-removability in Lemma 3(e), and that argument handles A = B separately and correctly.

Lemma 5 (mean square) is correct: the cross terms are geometric sums bounded by 2/|1−ω_kω̄_l|. The four ω are
distinct iff θ₀ ≠ π/2, since −e^{iθ₀} = e^{−iθ₀} iff θ₀ = π/2. Σ|a_k|² = r²/2. For A = B the merged coefficients
are A e^{±iπ/4}, so β_n = 2A cos(π/4 + nπ/2) and |β_n| = √2·A = r for all n. Lemma 6 is a second, quantitative
lower bound, not needed for Theorem A.

Side remark (not a gap): §5 says that for special θ₀, β_n "can vanish on arithmetic progressions". This is
asserted without an example or proof. It concerns only lim vs limsup, which the paper does not claim. Either mark it
as a remark or drop it.

## 4. Transfer theorem (F–S Thm VI.5): HOLDS

- Analytic in |z| < ρ: yes.
- Finitely many singularities on |z| = ρ: yes.
- Δ-analyticity: D = ⋂ z_jΔ(φ,R) with R < ρ₂/ρ. Each z_jΔ₀ removes its own outward ray. The intersection removes all
  four rays and lies in D(0,ρ₂), so D ⊂ Ω.
- Local form: valid on D ∩ D(z_j,δ). This set lies inside the connected set D(z_j,δ) minus the ray, where Lemma 3(d)
  holds by the identity theorem. δ < 2|L| keeps h_{3−i} analytic, and the ray from z₃ is at distance ≥ 2|L| from z₁.
- Scale functions: σ_j = c_j(s^{1/2} + (z_j/2)s^{3/2}) + polynomial, τ = s^{5/2}.
- Coefficients: Γ(−½) = −2√π, Γ(−3/2) = 4√π/3, and the correction 1 + α(α+1)/(2n) = 1 + 3/(8n) for α = ½. These
  give the bracket (1 + 3/(8n) − 3z_j/(4n)) and the error O(ρ^{−n}n^{−7/2}).
- Real form: the referee's script asserts that −(2√π)⁻¹Σc_jz_j^{−n} equals −√(2ρ)/(r√π)ρ^{−n}β_n to 30 digits for
  n ≤ 40 in all three cases.

## 5. Independent numerics: HOLDS (agreement)

Method: I_n is computed straight from the t-integral, not from Lemma 2. At each node the integrand
√(A²e^{2ηlog t}+B²e^{2ηlog(1−t)}) is expanded as a power series in η, and the coefficients are integrated with a
hand-written tanh–sinh rule at 50 digits. Halving the step changes nothing above 10⁻⁴⁸. H_n comes from the
referee's own series square roots.

| (A,B) | 1/ρ | \|I_n/n!\|^{1/n}, n = 20 / 40 | N_40 = ρ⁴⁰40^{3/2}I_40/40! vs one-term law | max_{8≤n≤40} \|I_n/n! − two-term law\|·ρⁿn^{7/2} | mean N_n², 20≤n≤40, vs Lemma 6 limit |
|---|---|---|---|---|---|
| (1,1) | 0.63662 | 0.50718 / 0.55406 | −0.97732 vs −1.00000 | 6.8 | 1.018 vs 2ρ/π = 1.000 |
| (1,3) | 0.52169 | 0.41898 / 0.45389 | −0.96534 vs −0.93514 | 7.8 | 0.685 vs ρ/π = 0.610 |
| (2,5) | 0.54990 | 0.40798 / 0.47663 | 0.82981 vs 0.80361 | 6.5 | 0.575 vs ρ/π = 0.579 |

How to read the table:
- The root test creeps toward 1/ρ, as the n^{−3/2} factor predicts.
- Amplitude and phase follow Theorem B, including near-cancellations. For (1,3) at n = 28, N = −0.030 against a
  law of 0.035; the two-term residual stays bounded there (6.6).
- The two-term residual ρⁿn^{7/2}|·| is bounded (≈ 2–8 for n ≥ 8). That is the predicted O(n^{−7/2}).
- The mean squares are near Lemma 6's limits; a 21-term Cesàro window plus the O(1/n) correction accounts for
  the 5–12% gaps.
- For (1,1) the referee's values agree with the author's `numerics_output.txt` to all printed digits, e.g.
  N_40 = −0.9773236968.
- The paper's next-term claim for the square, f_n/law = 1 + (3 − 3π(−1)ⁿ)/(8n) + O(n^{−2}), is consistent:
  n(f_n/law − 1) = 1.5593 vs 1.5531 at n = 39, and −0.907 vs −0.803 at n = 40, with the gap shrinking like 1/n.

## 6. The two "imprecisions in the paper's sketch": real, and the replacement HOLDS WITH FIX

(a) Real, but proof.md overstates it slightly. For j ≥ 1 the k_j are the Taylor coefficients of K − 1 as well as
of K, so the sketch's sentence is not false for j ≥ 1. The actual gap is the j = 0 coefficient of K − 1, namely
k₀ − 1 = (1+e^{−u})^{−2} − 1. The stated bound e^{−0·u}/0! = 1 is true for it but useless: it gives |n!(η_D(n)−1)| ≤ n!
with no gain. The needed bound is 0 ≤ 1 − k₀(u) ≤ 2e^{−u}, which follows from (1+x)^{−2} ≥ 1 − 2x. That is exactly
the missing ingredient.

(b) Real. "Cannot cancel the branch points" needs the reason, and the reason is R(𝓔) ≥ 2ρ > ρ.

Proposed sentence: correct in substance. One fix: "radius exactly ρ because each z_j is a simple zero of H's
radicand" gives only the lower bound on the growth (non-removability). It omits that there is no singularity
nearer than ρ. Also, H has two radicands, not one. Fixed sentence:

> Splitting $K=1+(K-1)$, the first part has Borel transform exactly $H$, whose Taylor series has radius
> exactly $\rho$: the radicands $A^2+B^2e^{-2z}$ and $B^2+A^2e^{-2z}$ have no zeros in $|z|<\rho$, and each
> $z_j$ is a simple zero of one of them, so it is a genuine branch point of $H$. The bounds
> $|k_j(u)|\le e^{-ju}/j!$ $(j\ge1)$ and $0\le1-k_0(u)\le2e^{-u}$ make the Borel transform of the rest
> analytic in $|\zeta|<2\rho$, so the Cauchy--Hadamard formula gives the limsup, whatever the interference
> of the four branch points.

---

**Verdict: yes. The paper's Theorem "Radius zero" (limsup |I_n/n!|^{1/n} = 1/ρ, and radius 0 in η and in d) is
proved for every rhombus A, B > 0, including A = B. Only the sketch's wording needs the fix in 6.**

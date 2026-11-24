# GRW2 — Radius zero at the diamond end, for every rhombus

**Target.** Paper `paper/ieee/main.tex`, Theorem "Radius zero" (`thm:zero`, line 1229). For the rhombus
|x/A|ᵖ + |y/B|ᵖ = 1 (A, B > 0) with perimeter P, put

  r = √(A²+B²),  η = 1/p − 1,  I(η) = P/(4r(1+η)),  I_n = I⁽ⁿ⁾(0)/n!,
  L = log(B/A),  ρ = √(L² + π²/4).

**Claim.** limsup_{n→∞} |I_n/n!|^{1/n} = 1/ρ; hence the Taylor series of P at p = 1 (in η or in d = p − 1)
has radius 0.

**Result.** Proved, for every A, B > 0, by hand (no computer-assisted step). The audit gap in
`R2_rhombus_ab/result.md` (that H keeps a genuine singularity at distance ρ, and that nothing cancels it)
is closed in Lemma 3 and Theorem A; the limsup needs no Darboux transfer and is immune to cancellation
between the conjugate points (Cauchy–Hadamard is exact). The leading asymptotic form, with amplitude and
oscillation, is then proved in Theorem B by the Flajolet–Odlyzko transfer, and the non-cancellation of the
four contributions is made quantitative in Lemma 6.

Labels: **PROVED** = every step written here. **NUMERICAL** = `numerics.py` (mpmath), output in
`numerics_output.txt`.

---

## 0. Notation

  R₁(z) = A² + B²e^{−2z},  R₂(z) = B² + A²e^{−2z},  h_i = √R_i / r,  H = h₁ + h₂,
  K(η,u) = (1+e^{−u})^{−2−η},  k_j(u) = (1+e^{−u})^{−2}(−log(1+e^{−u}))^j / j!.

The four points
  z₁ = L + iπ/2,  z₂ = L − iπ/2 = z̄₁,  z₃ = −L + iπ/2 = −z̄₁,  z₄ = −L − iπ/2 = −z₁
all have modulus ρ. Write z₁ = ρe^{iθ₀}, θ₀ = arg(L + iπ/2) ∈ (0, π) (θ₀ = π/2 iff A = B).
Put ρ₂ = √(L² + 9π²/4) > ρ.

Relation to the paper's notation: with q = (A²−B²)/(A²+B²), g(z) = √(cosh z + q sinh z) satisfies
r²(cosh z + q sinh z) = A²e^z + B²e^{−z}, so e^{−z/2}g(z) = √R₁(z)/r and e^{−z/2}g(−z) = √R₂(z)/r. Hence
the paper's H(z) = e^{−z/2}[g(z)+g(−z)] is exactly h₁ + h₂ (equal on ℝ, where all square roots are
positive, hence everywhere by continuation).

---

## 1. The integral representation

**Lemma 1 (Laplace form). PROVED.** For −1 < η < 1,
  I(η) = ∫₀^∞ e^{−u} H(ηu) K(η,u) du.

*Proof.* The first-quadrant arc is x = A t^{1/p}, y = B(1−t)^{1/p}, t ∈ [0,1]; since 1/p − 1 = η its speed
is (1+η)√(A²t^{2η} + B²(1−t)^{2η}), so P = 4(1+η)∫₀¹√(A²t^{2η}+B²(1−t)^{2η}) dt. Put t/(1−t) = e^l. Then
t^{2η} = [t(1−t)]^η e^{ηl}, (1−t)^{2η} = [t(1−t)]^η e^{−ηl}, dt = t(1−t) dl, t(1−t) = w(l) := e^l/(1+e^l)².
So
  I(η) = ∫_ℝ w(l)^{1+η/2} √(A²e^{ηl}+B²e^{−ηl}) / r dl.
w is even; fold l ↦ −l on (−∞,0) and use, for u ≥ 0, w(u)^{1+η/2} = e^{−u}e^{−ηu/2}(1+e^{−u})^{−2−η} and
e^{−ηu/2}[√(A²e^{ηu}+B²e^{−ηu}) + √(A²e^{−ηu}+B²e^{ηu})]/r = h₁(ηu) + h₂(ηu) = H(ηu). Convergence:
on ℝ, 0 < H(x) ≤ (A+B)(1+e^{−x})/r, so the integrand is ≤ C e^{−u}(1+e^{|η|u}), integrable for |η| < 1. ∎

**Lemma 1′ (K-expansion). PROVED.** K(η,u) = Σ_{j≥0} k_j(u)η^j (entire in η), and for u ≥ 0:
|k_j(u)| ≤ e^{−ju}/j! (as 0 < log(1+e^{−u}) ≤ e^{−u}); 0 ≤ 1 − k₀(u) = 1 − (1+e^{−u})^{−2} ≤ 2e^{−u}. Also
|∂_η^j K(η,u)| = K·log^j(1+e^{−u}) ≤ 1 for η ≥ −2, u ≥ 0. ∎

---

## 2. Smoothness and the exact coefficient formula

**Lemma 2 (coefficient formula). PROVED.** I ∈ C^∞(−1,1), and for every n ≥ 0
  I_n = Σ_{j=0}^{n} H_{n−j} M_{n,j},   H_m := [z^m]H(z),   M_{n,j} := ∫₀^∞ e^{−u}u^{n−j}k_j(u) du.
Moreover M_{n,0} = n!·η_D(n) with η_D(n) = Σ_{m≥0}(−1)^m(m+1)^{−n} (n ≥ 1; M_{0,0} = 1/2), and
  |M_{n,j}| ≤ (n−j)! / ( j!·(j+1)^{n−j+1} ).

*Proof.* Strip bound: for |Im z| < π/2, e^{−2z} has argument in (−π, π), so R_i(z) ∉ (−∞, 0] and the
principal roots are analytic in the strip S = {|Im z| < π/2}; |R₁(z)| ≤ A² + B²e^{−2Re z} ≤ (A + Be^{−Re z})²,
so |H(z)| ≤ (A+B)(1+e^{−Re z})/r. Cauchy's estimate on the disk |z − x| ≤ π/4 ⊂ S (x real) gives
  |H^{(m)}(x)| ≤ m!(4/π)^m · C₀ max(1, e^{−x}),  C₀ = 2e^{π/4}(A+B)/r.
Fix δ ∈ (0,1). For |η| ≤ δ, F(η,u) := H(ηu)K(η,u) satisfies, by Leibniz and Lemma 1′,
  |∂_ηⁿF| ≤ Σ_j C(n,j) u^{n−j} |H^{(n−j)}(ηu)| ≤ C_n (1+u)^n e^{δu},
and e^{−u}(1+u)^n e^{δu} is integrable. Dominated differentiation under the integral, n times, gives
I ∈ C^∞(−δ, δ) for each δ < 1 and I⁽ⁿ⁾(0) = ∫₀^∞ e^{−u}∂_ηⁿF(0,u) du, with
∂_ηⁿF(0,u) = Σ_j C(n,j) u^{n−j}H^{(n−j)}(0)·j!k_j(u) = n! Σ_j H_{n−j}u^{n−j}k_j(u). Divide by n!.
M_{n,0}: integrate (1+x)^{−2} = Σ_m (−1)^m(m+1)x^m, x = e^{−u}, termwise (for n ≥ 2 absolutely; for n = 0, 1
directly: 1/2 and log 2). The bound: Lemma 1′ and ∫₀^∞ u^{n−j}e^{−(j+1)u}du = (n−j)!/(j+1)^{n−j+1}. ∎

(This is the step at which no uniformity in n is needed: the identity holds for each n separately.)

---

## 3. The singularities of H: location, type, non-removability

**Lemma 3. PROVED.**
(a) The zeros of R₁ are exactly L + iπ(k+½), k ∈ ℤ; those of R₂ are −L + iπ(k+½). All are simple.
(b) The zeros of modulus < ρ₂ are exactly z₁, z₂ (for R₁) and z₃, z₄ (for R₂), all of modulus ρ.
(c) H continues analytically to the star-shaped (hence simply connected) slit disk
  Ω = D(0, ρ₂) ∖ ⋃_{j=1}^{4} {t z_j : t ≥ 1};
in particular H is analytic in D(0, ρ), so |H_m| ≤ C_σ σ^{−m} for every σ < ρ, with C_σ = max_{|z|=σ}|H|.
(d) (local form) Let s = 1 − z/z₁. On D(z₁, δ) ∖ {t z₁ : t ≥ 1}, for small δ > 0,
  h₁(z) = c₁ s^{1/2} ψ₁(z),  c₁ = (A/r)·(−2z₁)^{1/2},  ψ₁ analytic on D(z₁,δ), ψ₁ = 1 + (z₁/2)s + O(s²),
with principal branches of both square roots. By the symmetries H(z̄) = conj H(z) and A ↔ B (which
swaps h₁ ↔ h₂ and L ↦ −L), the same holds at every z_j with
  c₁ = (A/r)(−2z₁)^{1/2},  c₂ = c̄₁,  c₃ = (B/r)(−2z₃)^{1/2} = (B/r)(2z̄₁)^{1/2},  c₄ = c̄₃,
and ψ_j = 1 + (z_j/2)s_j + O(s_j²), s_j = 1 − z/z_j. If A ≠ B the other summand h_{3−i} is analytic at
z_j; if A = B then z₃ = z₁, z₄ = z₂, h₁ = h₂ and the two contributions add (H = 2h₁).
(e) (non-removability) Each z_j is a genuine singularity of H: H has no analytic extension to any
neighbourhood of z_j. Hence the Taylor series of H at 0 has radius of convergence exactly ρ.

*Proof.* (a) R₁(z) = 0 ⟺ e^{−2z} = −A²/B² ⟺ z = log(B/A) + iπ(k+½); R₁′ = −2B²e^{−2z} ≠ 0. Same for R₂.
(b) |±L + iπ(k+½)|² = L² + π²(k+½)², minimal (= ρ²) exactly for k ∈ {0, −1}, next value ρ₂².
(c) Ω is star-shaped about 0: if a segment [0,z] met a slit at t z_j (t ≥ 1), then z = s z_j with s ≥ t,
so z itself would lie on the slit. R₁, R₂ have no zeros in Ω (by (b) their only zeros in D(0,ρ₂) are the
slit endpoints). On a simply connected domain a zero-free analytic function has exactly two analytic
square roots; take the one positive at 0. On Ω ∩ S (convex, connected) it coincides with the principal
root (both are continuous roots agreeing at 0). So h₁, h₂, H continue to Ω, and agree with the functions of
Lemmas 1–2 on ℝ. Cauchy's estimate on |z| = σ gives the coefficient bound.
(d) Since B²e^{−2z₁} = −A², R₁(z₁ + w) = A²(1 − e^{−2w}) = A²·2w·φ(w), φ(w) = (1−e^{−2w})/(2w), analytic
and zero-free for |w| < π (the zeros of 1 − e^{−2w} are iπℤ), φ(0) = 1, φ(w) = 1 − w + O(w²). With
w = z − z₁ = −z₁ s: R₁ = A²(−2z₁)·s·φ(z−z₁). On the radial segment z = t z₁, t ↑ 1 (which lies in S since
Im z = tπ/2 < π/2), h₁ is the principal root; there s = 1 − t > 0 and arg R₁ → arg(−z₁) ∈ (−π, π) (as
Im(−z₁) = −π/2 ≠ 0), so for t near 1 the principal root factorizes:
h₁ = (A/r)(−2z₁)^{1/2}·(1−t)^{1/2}·φ^{1/2}. Put ψ₁ = φ(z−z₁)^{1/2} (principal; φ is near 1 on a small disk).
Both sides of (d) are analytic on the connected set D(z₁,δ) ∖ {t z₁ : t ≥ 1} and agree on the segment, so
they agree there. ψ₁ = 1 − w/2 + O(w²) = 1 + (z₁/2)s + O(s²). The other points follow by the stated
symmetries. When A ≠ B, h₂ is analytic on D(z₁,δ) for δ < 2|L| = |z₁ − z₃| (its singular points nearest
to z₁ are z₃ and z₄, at distances 2|L| and ≥ 2|L|).
(e) Suppose H extended analytically to a disk D around z_j. If A ≠ B, then h_i (the summand singular at
z_j) = H − h_{3−i} would extend analytically to D (after shrinking D so h_{3−i} is analytic there), and the
identity h_i² = R_i/r² would persist on D; but R_i has a simple zero at z_j, and the square of an analytic
function has only zeros of even order. If A = B, apply the same argument to h₁ = H/2. So z_j is singular.
The Taylor series of H at 0 converges on D(0,ρ) by (c); if its radius exceeded ρ, its sum would be an
analytic extension of H to a neighbourhood of z₁. ∎

**Type.** Each z_j is a square-root branch point, c_j ≠ 0, |c_j| = √(2ρ)·A/r (j = 1,2) or √(2ρ)·B/r (j = 3,4);
when A = B the merged coefficient is 2·√(2ρ)A/r = √(2ρ)·(A+B)/r.

---

## 4. The uniformity step: the K-part is exponentially smaller

**Lemma 4 (remainder bound, all n). PROVED.** Let E_n := I_n/n! − H_n. For every σ ∈ (0, ρ) and n ≥ 1,
  |E_n| ≤ C_σ (2σ)^{−n} [ 1 + ½ Σ_{j≥1} (2σ)^j/(j!)² ]  =: C′_σ (2σ)^{−n}.
Consequently the formal Borel transform 𝓑(ζ) := Σ_n I_n ζⁿ/n! equals H(ζ) + 𝓔(ζ) with
𝓔(ζ) = Σ E_n ζⁿ analytic in |ζ| < 2ρ.

*Proof.* By Lemma 2, E_n = H_n(η_D(n) − 1) + Σ_{j=1}^n H_{n−j}M_{n,j}/n!.
First term: η_D(n) − 1 = Σ_{m≥1}(−1)^m(m+1)^{−n} is alternating with decreasing terms, so
|η_D(n) − 1| ≤ 2^{−n}, and |H_n| ≤ C_σσ^{−n}; together ≤ C_σ(2σ)^{−n}.
Second term: n!/(n−j)! ≥ j! gives |M_{n,j}|/n! ≤ 1/((j!)²(j+1)^{n−j+1}) ≤ 2^{−(n−j+1)}/(j!)², so
  Σ_{j=1}^n |H_{n−j}||M_{n,j}|/n! ≤ Σ_{j≥1} C_σσ^{−(n−j)}2^{−(n−j+1)}/(j!)² = C_σ(2σ)^{−n}·½Σ_{j≥1}(2σ)^j/(j!)².
The series converges (it is ≤ I₀(2√(2σ))). Hence limsup|E_n|^{1/n} ≤ 1/(2σ) for every σ < ρ, i.e. ≤ 1/(2ρ). ∎

The bound is uniform in n: one constant C′_σ for all n ≥ 1, depending only on (A, B, σ).

---

## 5. The theorem

**Theorem A (radius zero). PROVED.** For every A, B > 0:
(i) limsup_{n→∞} |I_n/n!|^{1/n} = 1/ρ.
(ii) The Borel transform 𝓑(ζ) = Σ I_nζⁿ/n! has radius of convergence exactly ρ; on |ζ| = ρ its
singularities are exactly the square-root branch points z₁, …, z₄ (two points ±iπ/2 if A = B), and
𝓑 − H is analytic in |ζ| < 2ρ.
(iii) The Taylor series of P at p = 1, in η or in d = p − 1, has radius 0; P is C^∞ but not analytic at p = 1.

*Proof.* (i)–(ii): 𝓑 = H + 𝓔 as formal series (Lemma 4). H has radius exactly ρ (Lemma 3(e)) and 𝓔
has radius ≥ 2ρ > ρ. Hence 𝓑 has radius ≥ ρ; if it were > ρ, then H = 𝓑 − 𝓔 would have radius > ρ,
contradicting Lemma 3(e). So the radius is ρ, and by the Cauchy–Hadamard formula
limsup|I_n/n!|^{1/n} = 1/ρ. The singular points of 𝓑 on |ζ| = ρ are those of H, by (c)–(e) of Lemma 3.
(iii) By (i), along a subsequence |I_n/n!|^{1/n} → 1/ρ > 0, and (n!)^{1/n} → ∞, so |I_n|^{1/n} → ∞: the
Taylor series Σ I_nηⁿ has radius 0. Since P = 4r(1+η)I and P, I are C^∞ (Lemma 2), their Taylor series
are related by the formal identity Ser(I) = Ser(P)/(4r) · Σ_k(−η)^k; if Ser(P) had positive radius, so would
the product of two convergent series, a contradiction. In d: p = 1 + d, η = −d/(1+d) is analytic with
analytic inverse d = −η/(1+η); the Taylor series of a composition of C^∞ maps is the formal composition of
their Taylor series (Faà di Bruno), so a convergent d-series for P would compose with the convergent
series of d(η) to a convergent η-series. Finally, analyticity at p = 1 would make these series converge. ∎

**What the cancellation question does and does not affect.** The four points z₁,…,z₄ (two of them complex
conjugate to the other two) contribute oscillating terms that can partly cancel for individual n (see the
`bracket` column of the numerics: e.g. β₅ = −0.050 for (A,B) = (1,2)). This cannot change the limsup: the
Cauchy–Hadamard formula is an identity, and the radius is fixed by the existence of one non-removable
singularity at distance ρ, which Lemma 3(e) supplies. Cancellation is relevant only to finer statements;
Theorem B and Lemma 6 show it is never total along the whole sequence. What is **not** proved for A ≠ B is
that lim (rather than limsup) exists: β_n can come close to 0 along subsequences, and for special θ₀ it can
vanish on arithmetic progressions (when θ₀/π is rational), where the leading term drops out. For A = B,
|β_n| = r for every n and the limit exists.

---

## 6. The leading asymptotic form (amplitude and oscillation)

**Named theorem used.** Flajolet–Odlyzko transfer with finitely many dominant singularities:
Flajolet & Sedgewick, *Analytic Combinatorics* (CUP 2009), Theorem VI.5 (with Theorems VI.1 and VI.3);
originally Flajolet & Odlyzko, SIAM J. Discrete Math. 3 (1990) 216–240. Statement used: let f be analytic in
|z| < ρ with finitely many singularities ζ₁,…,ζ_m on |z| = ρ; suppose there is a Δ-domain
Δ₀ = Δ(φ, R) = {|z| < R, z ≠ 1, |arg(z−1)| > φ} (R > 1, 0 < φ < π/2) such that f is analytic in the indented
disk D = ⋂_j ζ_jΔ₀, and that, as z → ζ_j in D, f(z) = σ_j(z/ζ_j) + O(τ(z/ζ_j)) with σ_j a finite combination of
(1−z)^α's and τ(z) = (1−z)^β. Then [zⁿ]f = Σ_j ζ_j^{−n}([zⁿ]σ_j(z) + O(n^{−β−1})).
Coefficients: [zⁿ](1−z)^α = n^{−α−1}/Γ(−α)·(1 + α(α+1)/(2n) + O(n^{−2})) (α ∉ ℕ), and 0 for α ∈ ℕ, n > α.

**Hypotheses verified for f = H.** Analytic in |z| < ρ: Lemma 3(c). Singularities on |z| = ρ: exactly
z₁,…,z₄ (Lemma 3(b),(e)). Indented domain: take any φ ∈ (0, π/2) and R ∈ (1, ρ₂/ρ). A point t z_j, t ≥ 1,
has arg(z − z_j) = arg z_j, so it is excluded from z_jΔ₀; hence D ⊂ Ω ∩ D(0, Rρ) ⊂ Ω, where H is analytic
(Lemma 3(c)). Local expansions: by Lemma 3(d), as z → z_j in D,
  H(z) = c_j s_j^{1/2} + c_j(z_j/2)s_j^{3/2} + P_j(s_j) + O(|s_j|^{5/2}),
where P_j is a polynomial of degree ≤ 2 (the Taylor part of the analytic summand at z_j); the local
identity of Lemma 3(d) holds on all of D near z_j because D ∩ D(z_j,δ) is connected and contains the
radial segment. Take σ_j = c_j(s^{1/2} + (z_j/2)s^{3/2}) + P_j(s), τ = s^{5/2}.

**Theorem B (asymptotic form). PROVED.** As n → ∞,
  I_n/n! = −(1/(2√π)) n^{−3/2} Σ_{j=1}^{4} c_j z_j^{−n} (1 + 3/(8n) − 3z_j/(4n)) + O(ρ^{−n}n^{−7/2}),
and the leading term equals, in real form,
  I_n/n! = −√(2ρ)/(r√π) · ρ^{−n} n^{−3/2} · β_n + O(ρ^{−n}n^{−5/2}),
  β_n = B(−1)ⁿ cos((n−½)θ₀) − A sin((n−½)θ₀).
(For A = B the j = 1,3 and j = 2,4 points coincide; the formulas hold as written, and give the paper's
square law f_n ∼ √2 cos(π(3−2n)/4) n!(2/π)ⁿ n^{−3/2} when A = B = 1.)

*Proof.* Transfer (hypotheses above) with Γ(−½) = −2√π, Γ(−3/2) = 4√π/3:
[zⁿ]c_j s_j^{1/2} = −c_j z_j^{−n} n^{−3/2}(1 + 3/(8n) + O(n^{−2}))/(2√π), and
[zⁿ]c_j(z_j/2)s_j^{3/2} = c_j(z_j/2)z_j^{−n}n^{−5/2}(3/(4√π))(1+O(1/n)) = −c_jz_j^{−n}n^{−3/2}(−3z_j/(4n))/(2√π)·(1+O(1/n)).
P_j contributes nothing for n > 2, τ contributes O(ρ^{−n}n^{−7/2}). So H_n has the displayed expansion,
and I_n/n! = H_n + E_n with |E_n| ≤ C′_σ(2σ)^{−n} (Lemma 4, any σ ∈ (ρ/2, ρ)) is absorbed in the error.
Real form: θ₀ ∈ (0,π), so −2z₁ = 2ρe^{i(θ₀−π)} and 2z̄₁ = 2ρe^{−iθ₀} have principal roots
√(2ρ)e^{i(θ₀−π)/2} and √(2ρ)e^{−iθ₀/2}; z₃^{−n} = (−z̄₁)^{−n} = (−1)ⁿρ^{−n}e^{inθ₀}. Then
c₁z₁^{−n} + c₂z₂^{−n} = 2Re[...] = −2(A/r)√(2ρ)ρ^{−n} sin((n−½)θ₀),
c₃z₃^{−n} + c₄z₄^{−n} = 2(B/r)√(2ρ)(−1)ⁿρ^{−n} cos((n−½)θ₀).
Multiply by −n^{−3/2}/(2√π). ∎ (NUMERICAL: the real form equals the four-point sum to 35 digits for all
n ≤ 60 in all four test rhombi — an `assert` in `numerics.py`.)

**Lemma 5 (mean-square of exponential sums). PROVED.** Let ω₁,…,ω_m be distinct points of the unit circle,
a_1,…,a_m ∈ ℂ, and S_n = Σ a_kω_kⁿ. Then (1/N)Σ_{n=1}^N |S_n|² → Σ|a_k|². Hence limsup|S_n| ≥ (Σ|a_k|²)^{1/2}.
*Proof.* |S_n|² = Σ_k|a_k|² + Σ_{k≠l}a_kā_l(ω_kω̄_l)ⁿ, and |Σ_{n=1}^N(ω_kω̄_l)ⁿ| ≤ 2/|1−ω_kω̄_l| for k ≠ l. ∎

**Lemma 6 (the four contributions never cancel along the whole sequence). PROVED.**
β_n = Σ_k a_kω_kⁿ with
  ω ∈ {−e^{iθ₀}, −e^{−iθ₀}, e^{iθ₀}, e^{−iθ₀}},  a = {(B/2)e^{−iθ₀/2}, (B/2)e^{iθ₀/2}, −(A/2i)e^{−iθ₀/2}, (A/2i)e^{iθ₀/2}}.
If A ≠ B (θ₀ ≠ π/2) the four ω are distinct and (1/N)Σβ_n² → (A²+B²)/2 = r²/2. If A = B the ω merge in
pairs into ±i with merged coefficients of modulus A each (they add, they do not cancel) and |β_n| = √2A = r
for every n. Consequently
  limsup_n ρⁿ n^{3/2} |I_n/n!| ≥ √(ρ/π) (A ≠ B),   = √(2ρ/π) (A = B), and
  (1/N)Σ_{n≤N} (ρⁿn^{3/2}I_n/n!)² → ρ/π (A ≠ B),  2ρ/π (A = B).
*Proof.* Expand cos and sin; for A = B compute the merged coefficients directly:
i: (B/2)e^{iπ/4} − (A/2i)e^{−iπ/4} = A e^{iπ/4}; −i: its conjugate. Then Lemma 5 and Theorem B (whose
error is o(ρ^{−n}n^{−3/2})); the constant is (√(2ρ)/(r√π))²·r²/2 = ρ/π, resp. ·r² = 2ρ/π. ∎

Lemma 6 gives a second, quantitative proof of the lower bound in Theorem A(i) and answers the
cancellation question: the conjugate pair (and the reflected pair) can cancel at individual n, never in mean
square.

---

## 7. Numerical evidence (NUMERICAL) — `numerics.py`, output `numerics_output.txt`

Method: H_m by power-series square roots of R₁, R₂ at 45 digits; I_n = Σ_j H_{n−j}M_{n,j} (Lemma 2) with
M_{n,0} = n!η_D(n) exact (mpmath `altzeta`) and M_{n,j}, j ≥ 1, by tanh–sinh quadrature at 45 digits;
n ≤ 60. **Independent check** of the coefficients (not using Lemma 2's integrals): the optimally truncated
series Σ I_nηⁿ is compared with I(η) from direct quadrature of the arc-length integral at η = ±0.04; the
error must be of the size of the smallest term (≈ e^{−ρ/|η|} scale).

Run: `python numerics.py 60` (6.5 min, mpmath 1.4.1). Full table every 5th n in `numerics_output.txt`.
Normalised quantity: N_n := ρⁿ n^{3/2} I_n/n!; Theorem B predicts N_n = −√(2ρ)/(r√π)·β_n + O(1/n).

| (A,B) | ρ | 1/ρ | \|I_n/n!\|^{1/n} at n = 20 / 40 / 60 | N_60 vs law1 | max_{10≤n≤60} \|I_n/n! − law2\|ρⁿn^{7/2} | mean N_n² (30≤n≤60) vs Lemma 6 | \|E_n\|(1.5ρ)ⁿ at n = 60 | truncation vs direct quadrature, η = +0.04 / −0.04 |
|---|---|---|---|---|---|---|---|---|
| (1,1) | 1.570796 | 0.636620 | 0.5072 / 0.5541 / 0.5745 | −0.98548 vs −1.00000 | 4.62 | 1.014 vs 2ρ/π = 1.000 | 7.0e-11 | 1.7e-21 / 9.4e-22 (smallest term 5e-19) |
| (1,2) | 1.716932 | 0.582434 | 0.4540 / 0.4987 / 0.5261 | −1.03828 vs −1.04175 | 5.92 | 0.578 vs ρ/π = 0.547 | 7.2e-11 | 3.0e-22 / 6.7e-22 (1e-21) |
| (1,5) | 2.248931 | 0.444656 | 0.3568 / 0.3865 / 0.3991 | 0.71289 vs 0.72658 | 10.51 | 0.742 vs ρ/π = 0.716 | 5.1e-11 | 2.6e-28 / 7.0e-28 (1e-26) |
| (2,3) | 1.622283 | 0.616415 | 0.4720 / 0.5323 / 0.5554 | 0.89758 vs 0.90137 | 4.90 | 0.472 vs ρ/π = 0.516 | 6.2e-11 | 3.6e-21 / 9.2e-22 (1e-20) |

Reading the table.
* The root test |I_n/n!|^{1/n} creeps up to 1/ρ slowly (the n^{−3/2} and β_n factors cost ~ (3/2)log n/n), as expected; it is not, by itself, convincing — the next columns are.
* N_n follows the one-term law with relative error O(1/n) (e.g. (1,2): −1.0383 vs −1.0418 at n = 60), and the two-term law leaves a residual bounded by ≤ 10.5·ρ^{−n}n^{−7/2} over 10 ≤ n ≤ 60, i.e. the predicted O(n^{−7/2}) with a bounded constant, in all four cases. The amplitude √(2ρ)/(r√π) and the oscillation β_n (phase θ₀ = 0.3677π, 0.2461π, 0.4196π for (1,2), (1,5), (2,3)) are thus confirmed, including the near-cancellations (β₃₅ = 0.27 for (1,2), β₄₅ = −0.17 for (2,3)), where N_n is correspondingly small and still matches the law.
* The mean square over a 31-term window is within 6–10% of ρ/π (A ≠ B), resp. within 1.4% of 2ρ/π (square), consistent with Lemma 6 (the Cesàro error decays like 1/N times 1/|1−ω_kω̄_l|).
* E_n(1.5ρ)ⁿ → 0 geometrically, as Lemma 4 predicts (E_n = O((2σ)^{−n}) for every σ < ρ).
* Independent check: the optimally truncated series reproduces the direct arc-length quadrature of I(±0.04) to an error below the smallest term (1e-21…1e-28), at truncation orders 36–59. This tests the computed I_n themselves, not only Lemma 2's algebra.

---

## 8. Status table

| # | Step | Status |
|---|---|---|
| 1 | Laplace form I(η) = ∫₀^∞ e^{−u}H(ηu)K(η,u)du, −1 < η < 1; H = h₁ + h₂ equals the paper's e^{−z/2}[g(z)+g(−z)] | PROVED |
| 2 | I ∈ C^∞(−1,1); exact formula I_n = Σ_j H_{n−j}M_{n,j}; M_{n,0} = n!η_D(n); M-bound | PROVED |
| 3 | H analytic in D(0,ρ) and in the slit disk Ω up to ρ₂; singularities on \|z\| = ρ are exactly ±L ± iπ/2, simple square-root branch points with explicit c_j ≠ 0, non-removable; radius of H exactly ρ | PROVED |
| 4 | Uniform remainder bound \|E_n\| ≤ C′_σ(2σ)^{−n}, all n, every σ < ρ; 𝓑 − H analytic in \|ζ\| < 2ρ | PROVED |
| 5 | limsup\|I_n/n!\|^{1/n} = 1/ρ (Cauchy–Hadamard; immune to cancellation) | PROVED |
| 6 | Radius 0 for P in η and in d; P is C^∞, not analytic at p = 1 | PROVED |
| 7 | Leading asymptotic I_n/n! = −√(2ρ)/(r√π)ρ^{−n}n^{−3/2}β_n(1 + …) with two-term correction, error O(ρ^{−n}n^{−7/2}) — Flajolet–Odlyzko Thm VI.5, hypotheses checked | PROVED (O-constants not made explicit, depend on A,B) |
| 8 | Non-cancellation: mean square of ρⁿn^{3/2}I_n/n! → ρ/π (A≠B), 2ρ/π (A=B); limsup ≥ √(ρ/π) | PROVED |
| 9 | lim (not limsup) of \|I_n/n!\|^{1/n} | PROVED for A = B; OPEN for A ≠ B (not claimed by the paper) |
| 10 | Numerical agreement for (1,1), (1,2), (1,5), (2,3), n ≤ 60 | NUMERICAL |

**Verdict.** The theorem as stated in the paper (thm:zero) is **PROVED for every rhombus A, B > 0**,
including its "in particular" clause, by hand and without computer assistance.

**Notes on the paper's sketch (main.tex lines 1233–1238).** Nothing in the statement is wrong. Two
imprecisions in the sketch: (1) "the bound |k_j(u)| ≤ e^{−ju}/j! on the Taylor coefficients of K−1" —
the k_j are the coefficients of K, not K−1, and the j = 0 piece needs the separate bound
1 − k₀(u) ≤ 2e^{−u} (it is what produces the factor 2^{−n} through η_D(n) − 1 = O(2^{−n})); (2) "the rest
cannot cancel the branch points" should say *why*: its Borel transform is analytic in |ζ| < 2ρ. A sketch
sentence the paper can use verbatim:

> Splitting $K=1+(K-1)$, the first part has Borel transform exactly $H$, whose Taylor series has radius
> exactly $\rho$ because each $z_j$ is a simple zero of $H$'s radicand; the bounds
> $|k_j(u)|\le e^{-ju}/j!$ $(j\ge1)$ and $0\le1-k_0(u)\le2e^{-u}$ make the Borel transform of the rest
> analytic in $|\zeta|<2\rho$, so the Cauchy--Hadamard formula gives the limsup, whatever the interference
> of the four branch points.

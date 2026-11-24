# verify_proof — line-by-line audit of the proof in P4c `result.md` (Ahmed's Open Problem 4)

Auditor key: **proof**. Scratch work: `<scratch>/p4c_proof/` (scripts `s1`–`s6`, outputs `*_out.txt`).
Default stance was "wrong or incomplete until shown otherwise".

## 0. Verdict

**The theorem is true. The proof as written in `result.md` §2 has gaps; I closed them below.**
I re-derived every formula: R, R₃, R₄ (exactly, in sympy), the decomposition, the Mellin transforms and the
T_k lemma. All of them are correct. Two gaps are real ones in the logic:

* the Δ₁ remainder is written only as "+⋯", with no estimate and no uniformity in σ (gap G1);
* the e² coefficient of C, which R₄ needs, is never derived (gap G3).

There are also several write-up errors (G4–G10). After my closures:

* **R(n) and R₃(n) are proved.** The error term is an absolute O((1+log n)⁴/n⁴), valid for all real n ≥ 2.
* **R₄(n) and the O(log⁵n/n⁵) remainder of §1 are also proved.** `result.md` §6 still calls R₄ "NUMERICAL". It can now be called PROVED.
* **R_ρ in §5 is proved.** `result.md` wrote no proof for it; I supply one in §4.

---

## 1. Statement ready for the paper

Let A(n)=∬_{[0,1]²}√(1+n²x^{2n−2}y^{2n}+n²x^{2n}y^{2n−2}) dx dy. This is the area of the graph of z=xⁿyⁿ.
Put ℓ=log n, λ=log 2, ζ₂=π²/6. Then for real n→∞:

$$A(n)=1+\frac{\sqrt2}{n}+\frac{R(n)}{n^2}+\frac{R_3(n)}{n^3}+\frac{R_4(n)}{n^4}+O\!\Big(\frac{\log^5 n}{n^5}\Big)$$

The implied constant is absolute. Truncating after R₃ leaves an error of O((1+log n)⁴/n⁴) for every n ≥ 2. The coefficients are:

* R(n) = −½ℓ² + βℓ + γ, with β = λ/2 − 1 and γ = −(1+√2) + λ/2 − π²/12 − λ²/8.
* R₃(n) = ⅓ℓ³ + b₃ℓ² + c₃ℓ + d₃, with b₃ = (1−λ)/2, c₃ = λ²/4 − λ/2 + 1 + ζ₂, and
  d₃ = −λ³/24 + λ²/8 − λ/2 − ζ₂λ/2 + ζ₂/2 − ζ(3)/2 + 1 + 3√2/2.
* R₄(n) = −⅛ℓ⁴ + (λ/4 − 1/12)ℓ³ − (5/8 + π²/8 + 3λ²/16 − λ/8)ℓ²
  + (−5/4 − π²/24 − λ²/16 + λ³/16 + 5λ/8 + π²λ/8 + 3ζ(3)/4)ℓ
  + (−3π⁴/160 − 5/4 − 3√2/4 − 5π²/48 − 3λζ(3)/8 − π²λ²/32 − 5λ²/32 − λ⁴/128 + λ³/96 + π²λ/48 + ζ(3)/8 + 5λ/8).

All constants lie in ℚ[log 2, π², ζ(3), √2]. Euler's constant cancels identically.
The proof gives a structural statement as well: with N=√2 n, L=log N and e=1/n,

$$A(n)-1=\frac{\Phi(n)}{1+\frac1{2n}}-\frac{C(n)}{n^2},$$

* Φ is log-free.
* C(n) = poly₁ + e(poly₁ − poly₂ − 1/(2√2)) + e²(¾poly₁ − 5/4·poly₂ + ½poly₃ − 3/(2√2)) + O(e³(L+1)⁵).
* poly_k is the explicit polynomial of degree k+1 in L defined in Lemma 3.

Numerically: β = −0.6534264097…, γ = −2.9501636322…, b₃ = 0.1534264097…, c₃ = 2.4184737300…, d₃ = 2.4722752293…
I recomputed all of them to 45 digits (`s1_out.txt`); they match `result.md`.

## 2. Complete proof (my re-derivation; every step checked)

Notation: e=1/n, N=√2n, L=log N (≥1 for n≥2), f(w)=1+w−√(1+w²), g=f/w, W₀(σ)=Ne^{−σ}, and
T_k(φ)(N)=∫₀^N log^k(N/w) φ(w)dw = ∫₀^∞σ^k φ(W₀)W₀ dσ.

**Lemma 1 (exact decomposition) — `result.md` (a): correct.**
Substitute x=e^{−s/n}, y=e^{−t/n}, σ=s+t, δ=s−t, and use δ↔−δ. Then
A−1 = n^{−2}∫₀^∞e^{−eσ}∫₀^σ[√(1+W²)−1]dδdσ with W=W₀K, K=e^{eσ/2}√(cosh eδ).
The identity √(1+W²)−1 = W − f(W) holds. Both pieces are absolutely integrable, because 0≤f≤min(1,W).
* The W-piece is ∬ n x^{n−1}y^{n−1}√(x²+y²) dxdy. Putting y=ux on {y<x} turns it into (2n/(2n+1))Φ(n), with Φ=∫₀¹u^{n−1}√(1+u²)du.
* The f-piece gives C(n)/n².

Checks: decomposition against a direct (x,y) quadrature at n=3, 6, 12 agrees to ≤3e−27. The diagonal identity holds to 1e−36 (`s2_out.txt`).

**Lemma 2 (Φ) — correct.** Repeated integration by parts gives an exact finite expansion:
Φ = Σ_{j<J}(−1)^jφ^{(j)}(1)/(n(n+1)…(n+j)) + remainder, where φ=√(1+u²). The remainder is at most sup|φ^{(J)}|/(n…(n+J)).
In powers of e this gives:
* Φ = √2e − (√2/2)e² + (3√2/4)e³ − (7√2/8)e⁴ + (9√2/16)e⁵ + …
* Φ/(1+e/2) = √2e − √2e² + (5√2/4)e³ − (3√2/2)e⁴ + O(e⁵).

The expansion contains no log terms (sympy, `s1_out.txt`).

**Lemma 3 (the "Tauberian lemma") — `result.md` (c): correct. It is actually an exact identity.**
For N ≥ 1 and k ≥ 0:
$$T_k(g)(N)=\mathrm{poly}_k(N)+(-1)^k\frac{k!}{2N}+(-1)^{k+1}\theta\,\frac{k!}{8\cdot 3^{k+1}N^3},\qquad\theta\in[0,1],$$
where poly_k(N) = L^{k+1}/(k+1) + Σ_{i=0}^k k!/(k−i)!·m̃_{i+1}L^{k−i}, and m̃_j=[s^j]E_g(s)/(1−s), with E_g(s)=Γ(½+s/2)Γ(1−s/2)/√π.

*Proof.* Split the integral at w=1 and write g=1/w+(g−1/w) on [1,N]. Expand log^k(N/w)=(L−log w)^k binomially. This gives poly_k exactly, minus the tail ∫_N^∞log^k(N/w)(g−1/w)dw, with
c_i = ∫₀¹(−log w)^i g + ∫₁^∞(−log w)^i(g−1/w) = i!·[s^i]c(s), where c(s)=M_g(1−s)−1/s.
* c(s) is analytic for |Re s|<1.
* The identity c(s)=M_g(1−s)−1/s holds for 0<s<1 because ∫₁^∞w^{−s−1}=1/s.
* The tail is exact: g−1/w+1/(2w²) = 1/(2w²) − 1/(w(√(1+w²)+w)) lies in [0, 1/(8w⁴)]. This gives the k!/(2N) term and the θ-term.

*Mellin transform:* integrate by parts twice. The boundary terms vanish for 0<Re z<1, and
M_g(z) = (1/(z(1−z)))∫₀^∞w^z(1+w²)^{−3/2}dw = Γ((1+z)/2)Γ(1−z/2)/(√π z(1−z)) = E_g(z)/(z(1−z)).
This is a convergent beta integral, so no analytic continuation is needed. This route is cleaner than the "continuation" argument in `result.md`, and the two agree. Since E_g(1−s)=E_g(s), we get c(s)=[E_g(s)/(1−s)−1]/s, which gives μ_i=m̃_{i+1}.

*Expansion of E_g:* log E_g(s) = −λs + ζ₂s²/2 − ζ(3)s³/4 + … The γ_E terms cancel between ψ(½) and ψ(1).
So m̃₁=1−λ, m̃₂=1−λ+λ²/2+ζ₂/2, and m̃₃=m̃₂−ζ(3)/4−λζ₂/2−λ³/6, exactly as claimed.

Checks:
* c₀…c₃ from direct quadrature equal i!m̃_{i+1} to ≤1e−38.
* M_g(0.5) and M_g(0.4+0.7i) agree with the formula to ≤2e−18.
* θ_k lies in [0,1] for k=0..4 at N=1, 2, 14.1, 1000, and θ→1 (`s2_out.txt`).

**Lemma 4 (pointwise majorants) — new, needed to close G1/G3.** Define ψ(w)=wf′(w)=w/(√(1+w²)(√(1+w²)+w)), χ₁=wψ′ and χ₂=wχ₁′. Then:
* 0 ≤ ψ(w) ≤ min(w, 1/(2w));
* |χ₁| ≤ ψ and |χ₂| ≤ ψ.

*Proof.* Put w=sinh t and T=tanh t. Then ψ=tanh t·e^{−t} and w d/dw = T(1−T²)d/dT, so
χ₁/ψ = 1−T−T² and χ₂/ψ = 1−3T(1−T)(1+T)². Both lie in [−1,1] on [0,1].
The bound ψ ≤ 1/(2w) is equivalent to 2sinh²t ≤ cosh²t + sinh t cosh t (`s3_out.txt`).

Two consequences follow. Let κ(σ,δ)=log K = eσ/2 + ½log cosh(eδ); then 0≤κ≤eσ for 0≤δ≤σ.
* For u∈[W₀,W₀e^κ]: ψ(u), |χ_j(u)| ≤ e^κ m(σ), where m(σ)=min(W₀,1/(2W₀)).
* e^{−eσ}e^κ ≤ 1. **This is the uniformity mechanism.** It also covers the range σ≳n, where K−1 is not small.

Put I_k=∫₀^∞σ^k m dσ. Then I_k ≤ L^k/2 + Σ_i k!/(k−i)!L^{k−i} ≤ (k!+½)(L+1)^k, and T_k(g)(N) ≤ (k!+1)(L+1)^{k+1}.

**Lemma 5 (Δ₁) — `result.md` (b)/(d): result correct, proof incomplete (G1, G3). Closed here.**
Δ₁ = ∫e^{−eσ}∫₀^σ[f(W₀e^κ)−f(W₀)]dδdσ. Taylor-expand in κ (not in K−1). With Ψ(τ)=f(W₀e^τ) we have Ψ′=ψ, Ψ″=χ₁, Ψ‴=χ₂, evaluated at W₀e^τ.

(i) Δ₁ = e·T₁(g)(N) + r₁ with |r₁| ≤ (13/12)e²I₃ = O(e²(L+1)³).
* The Taylor remainder is ≤ (e²/2)I₃.
* ∫₀^σκdδ = eσ²/2 + k₁ with 0≤k₁≤e²σ³/12, because log cosh x ≤ x²/2.
* |e^{−eσ}−1| ≤ eσ.
* Finally (e/2)∫σ²ψ(W₀)dσ = (e/2)T₂(f′) = eT₁(g), by the integration by parts T_{k+1}(f′)=(k+1)T_k(g). The boundary terms vanish since f(w)~w at 0.

(ii) Δ₁ = eT₁ + e²[¾T₁(g) − (5/4)T₂(g)] + r₂, with |r₂| ≤ ⅔e³I₄ + 0.0563e⁴I₅ + e⁵I₆/120 = O(e³(L+1)⁴).
The extra inputs are: x²/2−x⁴/12 ≤ log cosh x; ∫₀^σκ² = e²σ³/4 + O(e³σ⁴+e⁴σ⁵); and the identities
∫σ³ψ(W₀) = 3T₂(g), ∫σ³χ₁(W₀) = 6T₁(g).
Together these give e²∫σ³[(1/12−1/2)ψ + (1/8)χ₁] = e²[¾T₁ − (5/4)T₂].
**`result.md` never derives this e² coefficient** (gap G3). It agrees with its §2(d) bracket.

**Lemma 6 (Δ₂) — correct, truncation bound added (G2).** Since f ≥ 0 and the Taylor remainders of e^{−x} alternate in sign:
Δ₂ = −eT₂ + r₃ with 0 ≤ r₃ ≤ (e²/2)T₃; and Δ₂ = −eT₂ + (e²/2)T₃ + r₄ with −(e³/6)T₄ ≤ r₄ ≤ 0.
**The L⁵ in O(e³L⁵) comes only from this term.** The Δ₁ errors are one log smaller, because ψ is localized near σ=L, whereas g~1/w is not.

**Assembly — `result.md` (d)/(e): correct (exact sympy check).** Combine C=C₀+Δ₁+Δ₂ with C₀=T₁(g) and Lemmas 3, 5 and 6, and insert into E=Φ/(1+e/2)−e²C. The tails −e/(2√2) (from C₀) and −3e²/(2√2) (from eT₁ and −eT₂) enter at the right orders. Substituting L=ℓ+λ/2 reproduces β, γ, b₃, c₃, d₃ **exactly**, and also reproduces the R₄ formula of §2(d). All coefficient differences are 0 in `s1_out.txt`.
* Error after R₃: e²·O(e²(L+1)⁴) + O(e⁴).
* Error after R₄: e²·O(e³(L+1)⁵) + O(e⁵).

## 3. Gaps and errors found (with status)

| # | Where | Problem | Severity | Status |
|---|---|---|---|---|
| G1 | §2(b) Δ₁ | "f(W₀K)−f(W₀)=W₀f′(W₀)(K−1)+⋯". There is no remainder bound, and Taylor in K−1 is not uniform: K−1 is not small once σ≳n, and the region σ∈[0,∞) is unbounded. The O(e²L³) is asserted, not proved. | substantive (the R₃ proof depends on it) | **closed** (Lemmas 4, 5(i)) |
| G2 | §2(b) Δ₂ | "exact, termwise" is fine (dominated convergence: f≤W₀ and e<1), but no truncation bound is given. | minor | **closed** (Lemma 6) |
| G3 | §2(d) | The e² coefficient ¾poly₁ − 5/4·poly₂ of Δ₁ is stated without derivation, and the O(e³L⁵) error is unproved. | substantive for R₄ | **closed** (Lemma 5(ii), Lemma 6) |
| G4 | §1 vs §6 | §1 labels the expansion through R₄ + O(log⁵n/n⁵) "PROVED"; §6 lists R₄ as "NUMERICAL". | inconsistency | resolved: now PROVED |
| G5 | §3 | "c₃ = λ²/4 − λ/2 + 1 + 2B_D" is **false**: it gives 2.2058, not 2.4185. The correct prev-route relation is c₃ = (λ/2)² − λ(A_D−½) + 1 + 2B_D. The prev route's d₃ had a sign error in prev/step_11, 12, fixed in prev/step_13; the fixed version equals d₃ exactly. The closed forms of A_D, B_D, C_D and M_h are correct (sympy + quadrature). | write-up error | §3 is not needed for the §2 proof; **drop it from the paper or fix it** |
| G6 | §2(a) | "The W−1 part … gives ∬x^{n−1}y^{n−1}√(x²+y²)" omits a factor n. Also, the −1 belongs with f, since √(1+W²)−1 = W − f(W). The boxed formula is right. | slip | fix text |
| G7 | §2(c) | "Tauberian lemma" is a misnomer: it is an exact identity with an explicit O(N⁻³) term. The Mellin step "beta + analytic continuation" needs one sentence of justification; the IBP route above avoids it. | presentation | fix text |
| G8 | §5 | ρ=½ evidence "0.392, 0.086, 0.063, 0.062 (n=1e3,1e4,2e3→1.6e4)" misreports the run. The transcript values at n=2000, 4000, 8000, 16000 are 0.2538, 0.1613, 0.1008, 0.0621. | reporting error | fix text (decay is genuine; I reproduced it) |
| G9 | §5 | R_ρ is called "PROVED (same machinery)" but no proof is written. | gap | **closed** (§4 below) |
| G10 | §7 / citations | File names are offset: for N≥12, `step_N.py` carries the header "step_{N−2}.py", and `result.md` cites the header numbers. Examples: "step_23 comprehensive check" is file `step_25.py`; "step_10/11 Mellin" are files `step_12/13.py`. | reproducibility | rename or re-cite |
| G11 | §2(c) | poly_k uses L=log N, not log n. This is implicit. | notation | state explicitly |

I also checked the J_ρ series identity in §5. It is not needed for any asymptotic statement, and I did not audit it; only the author's numerical verification exists.

## 4. Closing G9: proof of the ρ-generalization (fixed ρ>0, m=ρn real)

Write E_ρ = 𝒥 − C_ρ with:
* 𝒥 = ∬x^{n−1}y^{m−1}√(n²y²+m²x²) = (1/(ρn))∬e^{−a−b}F(a/n, b/m) da db, where F(α,β)=√(e^{−2β}+ρ²e^{−2α});
* C_ρ = (1/(ρn²))∬_{s,t≥0}e^{−s/n−t/m}f(W) ds dt, with W = W₀K_ρ, W₀ = N_ρe^{−σ}, N_ρ = n√(1+ρ²);
* K_ρ² = (e^{2s/n}+ρ²e^{2t/m})/(1+ρ²) ∈ [1, e^{2max(s/n,t/m)}].

**𝒥 part.** Every second partial derivative of F is bounded by 3F ≤ 3√(1+ρ²). For example, ρ²e^{−2α}/F ≤ F. Hence
𝒥 = √(1+ρ²)/(ρn) + c₁/(ρn²) + O_ρ(n⁻³), with c₁ = −(1+ρ³)/(ρ√(1+ρ²)).

**C_ρ part.** Let x=s/n+t/m and κ=log K_ρ ≤ x. Then e^{−x}e^κ ≤ 1, so
|e^{−x}f(W)−f(W₀)| ≤ x·[m(σ)+min(1,W₀)], with x ≤ eσ·max(1,1/ρ).
Integrating over the measure σdσ gives ρn²C_ρ = T₁(g)(N_ρ) + O_ρ(e(L+1)³). Lemma 3 then yields R_ρ(n) = c₁/ρ − poly₁(N_ρ)/ρ, which is exactly the formula in `result.md`, with
E_ρ = √(1+ρ²)/(ρn) + R_ρ/n² + O_ρ(log³n/n³).

Numerical check (`s5_out.txt`). I computed D_ρ := ρn²C_ρ − T₁(g)(N_ρ) and the 𝒥-remainder independently:

| ρ | n | D_ρ·n/L³ | n²·(𝒥 remainder) | n²E_ρ − lead − R_ρ |
|---|---|---|---|---|
| ½ | 250 | −0.583 | 4.304 | 0.8717 |
| ½ | 1000 | −0.553 | 4.330 | 0.3919 |
| ½ | 4000 | −0.535 | 4.336 | 0.1613 |
| ½ | 16000 | −0.524 | 4.337 | 0.0621 |
| 2 | 1000 | −0.271 | 2.167 | 0.0635 |
| 2 | 16000 | −0.260 | 2.169 | 0.0094 |

Both remainders are bounded, as Lemma 3/4 predict, and the residual decays like log³n/n.

## 5. Independent checks performed (not a repeat of the planner's constant check)

1. **Exact assembly** (`s1_symbolic.py`). Sympy series of log E_g, Φ by exact IBP, and the C-expansion built from Lemmas 3, 5, 6. R, R₃, R₄ match `result.md` with zero coefficient differences. Euler's γ cancels.
2. **Decomposition against the raw definition** (`s2_lemma.py`). A direct 2-D (x,y) quadrature of A(n)−1 at n=3, 6, 12 agrees with Φ/(1+1/2n) − C/n² to ≤3e−27.
3. **Lemma 3** (`s2_lemma.py`). θ_k ∈ [0,1] for k≤4 and several N; c_i = i!m̃_{i+1} to 1e−38; Mellin M_g and M_h checked.
4. **Lemma 4** (`s3_bounds.py`). Exact rational forms χ_j/ψ, with sup|χ₁/ψ| = sup|χ₂/ψ| = 1.
5. **Remainder orders** (`s4_remainders.py`, 34 digits, n=10…10⁵):

   | n | r₁/(e²(L+1)³) | r₂/(e³(L+1)⁴) | (C − 3-term expansion)/(e³L⁵) |
   |---|---|---|---|
   | 10 | −0.290 | 0.135 | −0.155 |
   | 1000 | −0.324 | 0.119 | −0.0348 |
   | 100000 | −0.345 | 0.126 | −0.0303 |

   * r₁/(e²(L+1)³) moves slowly toward −5/12, the leading coefficient of X.
   * The last column tends to −1/30, the predicted leading term −L⁵/30 of C⁽³⁾.
   * (T_k quadrature − identity)·N³ → 1/72 and −1/108 for k=1, 2, which means θ → 1.
6. **Independent E(n)**. My own quadrature at n=10³, 10⁴, 10⁵ matches the `result.md` table to 5e−30, 4e−33 and 8e−35.
7. **R₄ closed form against data** (`s6_R4check.py`). r₄ = n⁴(E − √2/n − R/n² − R₃/n³) − R₄ takes the values 0.70, 0.44, 0.25, 0.071 at n=10³, 3·10³, 10⁴, 10⁵. Also r₄·n/ℓ⁵ → 1/30 (0.0351 at 10⁵), the leading coefficient of R₅.

## 6. Bottom line for the paper

You can quote the theorem in §1 above as proved, including R₄ and the O(log⁵n/n⁵) remainder, provided the paper's proof includes:
* Lemma 3 with the exact tail and the IBP Mellin evaluation;
* Lemma 4 (the majorants ψ ≤ min(w,1/2w) and |χ_j| ≤ ψ);
* Lemma 5 (the Taylor expansion in κ=log K, together with e^{−eσ}e^κ ≤ 1);
* Lemma 6.

Do not use the `result.md` text verbatim: it lacks the remainder estimates (G1, G3) and contains the false §3 relation (G5). The ρ-generalization can be stated as a proved proposition with error O_ρ(log³n/n³), using the argument in §4.

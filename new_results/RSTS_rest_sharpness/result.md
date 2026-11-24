# RSTS: the rest exponent √(2π) in obs:Mk is sharp

**Verdict: SHARP PROVED.** The rest has an explicit leading term of exact size N^{-1}e^{-√(2πN)}, with nonzero amplitude √π. The proof is by hand. The error term O(N^{-1/2}) in the leading law is proved but its constant is not computed (the theorem is not effective in that term). The other part of the rest (U_k below) has an explicit, arb-certified constant.

Two side findings:
1. The phase of the q=2 term is −π/4, not 0. RST1/result.md says "rest ~ (√π/N)e^{-√(2πN)}cos√(2πN)". That is wrong. The correct law is cos(√(2πN) − π/4). D13c had the phase π/4 right.
2. The proved bound 1030·N^{-1/2}e^{-√(2πN)} is sharp in the exponent but not in the power. The true size is N^{-1}.

## 0. What exactly is claimed

The display eq:Mkfull reads "first line with O(N^{-1}) inside the bracket, plus O(e^{-√(2πN)})". As written, the last term says nothing: the first line already carries an error N^{-2}e^{-√(πN)}, which is far larger than e^{-√(2πN)}. The same holds if the q=1 series is truncated at any finite order m: its remainder N^{-1-m/2}e^{-√(πN)} still swamps e^{-√(2πN)}.

So sharpness only makes sense for the exact split used in the proofs (CX6 §3, RST1). There M_k^main is an exact integral, not a truncated series:

  M_k^main = (2/π) Im ∫_0^∞ H(s) e^{iπs/2} s^{-1} (s/(1+s))^N ds,  H(s) = ½B(s/2, 3/2),  N = k+1,

and rest_k := M_k − M_k^main. The paper's text says this in words ("the main contribution of the branch points ±iπ/2"). The open problem is then: is √(2π) the right exponent for rest_k? The answer is yes, in the strongest form (a).

## 1. Theorem

Put N = k+1, θ₂ = √(2πN), λ = √(N/π), Λ = √(πN).

**Theorem (PROVED).** rest_k = S_k + U_k, where

  S_k = (2/π) Im ∫_0^{∞} H(s) e^{iπs/2} A^{[2]}_k(s) ds,  A^{[2]}_k(s) = −2(s+1)s^{k−1}/(s+2)^{k+1}

(path: the ray arg s = π/4 from 0, or the positive real axis; they agree), and

(i) as N → ∞,

  S_k = (√π/N) e^{−√(2πN)} [ cos(√(2πN) − π/4) + O(N^{−1/2}) ];

(ii) for all N ≥ 100,

  |U_k| ≤ C_U e^{−(3/2)√(π(N−3))},  C_U = 10 583 366.83 (arb upper bound; s4_constant_U.log).

Since (3/2)√π = √(9π/4) = 2.6587 > √(2π) = 2.5066, U_k is exponentially smaller than S_k.

**Constants of the q=2 term.** A₂ = √π = 1.772453850905516027298167483341145182798 (nonzero), β = −1, φ₂ = −π/4. In complex form, rest_k = (√π/N)e^{−θ₂} Re[e^{i(θ₂−π/4)}(1+O(N^{−1/2}))].

**Corollary (sharpness, PROVED).**

  limsup_{N→∞} N e^{√(2πN)} |M_k − M_k^main| = √π,  liminf_{N→∞} N e^{√(2πN)} |M_k − M_k^main| = 0.

So no bound C·N^b·e^{−c√N} with c > √(2π) holds for any C, b. The exponent √(2π) is sharp, and so is the power N^{−1}.

**First correction (DERIVED, not proved; numerically confirmed).** S_k = (√π/N)e^{−θ₂} Re[e^{i(θ₂−π/4)}(1 + ϖ₂/√N + O(N^{−1}))], with

  ϖ₂ = √2[(12 + 33π + 2π²) + i(12 − 33π + 2π²)]/(48√π) = 2.250894952878746262386866788 − 1.195718924738879428434810478 i.

The unjustified step is only that the formal second-order saddle expansion (s3_varpi2.py) is not backed by a written error bound. The same code applied to the q=1 integral gives the paper's ϖ₁ exactly (s3b_check_varpi1.log, difference 0 at 167 digits).

## 2. Proof

### 2.1 Exact decomposition (from CX6, PROVED there)

CX6 (8) gives, with Γ: s = ½ + rω, ω = e^{iπ/4}, K = 1/(2cos(πs/2)), E = e^{iπs/2}, A_k = P_k − Q_k, Q_k = s^k/(1+s)^{k+1}:

  rest_k = −J_k + (2/π) Im{ ∫_Γ H K A_k + ∫_Γ H(K−E)Q_k − ∫_0^{1/2} H E Q_k }.

Split A_k = A^{[2]}_k + A^{[≥3]}_k and K = E + (K−E) on the A^{[2]} part. Move ∫_Γ H E A^{[2]} to the ray from 0. The integrand is analytic in Re s > 0 and is O(s^{k−2}) at 0. The strip between Γ and the ray meets the real axis only in [0, ½], so ∫_Γ = ∫_ray − ∫_0^{1/2}. Hence rest_k = S_k + U_k with

  U_k = −J_k + (2/π) Im{ I₁ + I₂ + I₃ − ∫_0^{1/2} H E Q_k − ∫_0^{1/2} H E A^{[2]}_k },
  I₁ = ∫_Γ H(K−E)A^{[2]}_k,  I₂ = ∫_Γ H K A^{[≥3]}_k,  I₃ = ∫_Γ H(K−E)Q_k.

**Where A^{[2]} comes from.** CX6 (10): A_k(s) = (1/k!)∫_0^∞ e^{−(s+1)L}[(sL−β)^k/(1+w) − (sL)^k] dL, with w = e^{−L} and β = log(1+w). Expand (sL−β)^k binomially. The i=0 term is (sL)^k(1/(1+w) − 1) = (sL)^k(−w + w²/(1+w)). The i=1 term is −k(sL)^{k−1}β/(1+w) = −k(sL)^{k−1}(w + r₁), where r₁ := β/(1+w) − w. Keep only the "−w" parts: (1/k!)∫ e^{−(s+2)L}[−(sL)^k − k(sL)^{k−1}] dL = −s^k/(s+2)^{k+1} − s^{k−1}/(s+2)^k = A^{[2]}_k(s). This is also the m=1 term of the binomial series of (1+t)^{−1−ε} in CX6's Φ. It is the pole s = −2, i.e. the second image of ±iπ/2.

### 2.2 Bound on A^{[≥3]} (PROVED, by hand)

For |arg L| ≤ π/4, CX6 (11) gives |1+vw| ≥ ½ (0 ≤ v ≤ 1) and |β| ≤ 2|w|. Also β − w = −w²∫_0^1 t/(1+tw) dt, so |β − w| ≤ |w|². Hence |w²/(1+w)| ≤ 2|w|² and |r₁| = |β − w − w²|/|1+w| ≤ 4|w|². Extract e^{−2L} or e^{−iL} analytically. Then rotate each L-integral to arg L = −arg(s+3) or −arg(s+i+1), exactly as in CX6 §5 (same wedge, no singularities crossed, arcs vanish). This gives

  |A^{[≥3]}_k(s)| ≤ 2ρ^k/|s+3|^{k+1} + 4ρ^{k−1}/|s+3|^k + 2 Σ_{i=2}^k (2^i/i!) ρ^{k−i}/|s+i+1|^{k−i+1},  ρ = |s|.   (B3)

### 2.3 Bound on U_k (PROVED by hand; constants CERTIFIED with arb)

Use CX6 (13) on Γ: |H| ≤ 2, |K| ≤ 2e^{−ar}, |K−E| ≤ 2e^{−3ar}, a = π/(2√2). Use CX6 (14): F(j,q;λ) ≤ (4/λ)e^{λ/2+λq/√2} e^{−√(3λjq/√2)}.

- I₁: |A^{[2]}| ≤ 2ρ^{k−1}/|s+2|^k, since |s+1| ≤ |s+2|. So |I₁| ≤ 8F(k−1, 2; 3a), with decay e^{−√(9π(N−2)/2)}.
- I₂: by (B3), |I₂| ≤ 4[2F(k,3;a) + 4F(k−1,3;a) + 2Σ_{i≥2}(2^i/i!)F(k−i, i+1; a)]. Here 3a/√2 = 3π/4. The first two terms decay like e^{−√(9π(N−2)/4)}. For 2 ≤ i ≤ k/2, (k−i)(i+1) − 3(k−2) = (i−2)(k−3−i) ≥ 0, so each term is ≤ e^{−√(9π(N−3)/4)}, with coefficient sum ≤ e^{x}, x = 2e^{π/4}. For i > k/2, use CX6 (16): the sum is ≤ e^{ex}e^{−k/2}.
- I₃: CX6 (17): (2/π)|I₃| ≤ T₀ e^{−√(9πk/4)}.
- Elementary terms: |J_k| ≤ (log 2)^N/N! ≤ e^{−(N−1)/2}. The two [0,½] segments are ≤ (2/(πk))3^{−k} + (6/π)·2^{−2k}/(k−1), using H(s) ≤ 1/s on (0,½].
- Every exponent is ≥ (3/2)√(π(N−3)). For the e^{−k/2} and e^{−(N−1)/2} terms this holds for N ≥ 100: at N = 100 it is 49.5 ≥ 26.18, and the difference increases, since d/dN[(N−1)/2 − 1.5√(π(N−3))] > 0 once √(N−3) > 1.5√π.

Adding the arb-evaluated prefactors (s4_constant_U.py) gives C_U < 10 583 367. The constant is crude, dominated by e^{ex}. The true size is |U_k| ≈ (0.02–0.7)·(√π/N)e^{−√(3πN)} for 100 ≤ k ≤ 400 (table, column u).

### 2.4 The q=2 saddle (PROVED; O-constant not computed)

Write Z_k = ∫_ray H E A^{[2]}_k ds, so S_k = (2/π) Im Z_k. Put s = 2λrω, r > 0. Then e^{iπs/2} = e^{iΛrω}. Put x = 2/s = 1/(λrω). Then A^{[2]} = −(2/s)(1 + x/2)(1+x)^{−N}, and

  (1+x)^{−N} = e^{−Nx} e^{Nx²/2} e^{−Nτ(x)},  τ(x) = log(1+x) − x + x²/2,  |τ(x)| ≤ |x|³/(3(1−|x|)).

Here Nx = Λω̄/r and Nx²/2 = −iπ/(2r²) exactly. Stirling for the ratio Γ(z)/Γ(z+3/2), uniform in |arg z| ≤ π/4 (DLMF 5.11.13), gives H(s) = √(π/2)s^{−3/2}(1+η(s)) with |η(s)| ≤ c_H/|s| for |s| ≥ 1. Collect the factors:

  H E A^{[2]} ds = λ^{−3/2} p(r)(1+ε_N(r)) e^{Λψ(r)} dr,
  ψ(r) = iωr − ω̄/r,  p(r) = −(√π/2) ω^{−3/2} r^{−5/2} e^{−iπ/(2r²)},  1+ε_N = (1+η)(1+x/2)e^{−Nτ(x)}.

Facts about ψ:
- Re ψ(r) = −(r + 1/r)/√2 = −√2 − (r−1)²/(√2 r).
- ψ(1) = −√2 + i√2, so Λψ(1) = −θ₂ + iθ₂.
- ψ′(1) = 0 and ψ″(1) = −2ω̄, with Re(−ψ″(1)) = √2 > 0.
- p and ψ do not depend on N.

(a) **Region r ≥ ¼, λ ≥ 8.** Here |x| ≤ 4/λ ≤ ½ and |s| ≥ λ/2. Also N|τ| ≤ (2/3)π/(λr³) ≤ 134/λ. So |ε_N(r)| ≤ C/λ uniformly. The perturbation is therefore ≤ (C/λ)·λ^{−3/2}e^{−θ₂}∫|p|e^{−Λ(r−1)²/(√2r)}dr = O(N^{−3/2}e^{−θ₂}).

(b) **Laplace lemma for the fixed pair (p, ψ).**

  ∫_{1/4}^∞ p e^{Λψ} dr = e^{Λψ(1)} √(2π/(−Λψ″(1))) [p(1) + O(Λ^{−1})],

with the principal root. Proof: Off |r−1| ≤ Λ^{−2/5}, the modulus is ≤ e^{−√2Λ}e^{−cΛ^{1/5}} (and ≤ e^{−√2Λ}e^{−Λ(r−1)/(2√2)} for r ≥ 2), times the integrable |p|. Inside, put t = r−1 and expand ψ − ψ(1) = ψ₂t² + ψ₃t³ + O(t⁴) and p = p₀ + p₁t + O(t²). Since Λ|t|³ ≤ Λ^{−1/5}, e^{Λ(ψ₃t³+O(t⁴))} = 1 + Λψ₃t³ + O(Λt⁴ + Λ²t⁶). The odd terms p₁t and p₀Λψ₃t³ integrate to 0 against e^{Λψ₂t²}. The rest is O(Λ^{−3/2}), because Re ψ₂ = −1/√2. Finally ∫e^{Λψ₂t²}dt = √(π/(−Λψ₂)) up to e^{−cΛ^{1/5}}. (This is Olver, *Asymptotics and Special Functions*, Ch. 4 §6, in the form needed.)

(c) **Region 0 < r < ¼ (|s| < λ/2).** Here |1+2/s|² ≥ 1 + 2√2/|s| ≥ 1 + 4√2/λ. With log(1+u) ≥ u − u²/2 and N/λ² = π, this gives |s/(s+2)|^N ≤ e^{−2√2Λ+8π}. With |H| ≤ √2/|s| and |E| ≤ 1, the piece |s| ∈ [1, λ/2] is ≤ 2√2λe^{−2√2Λ+8π}. The piece |s| ≤ 1 is ≤ 4√2/(2^{k+1}(k−1)). Both are o(N^{−3/2}e^{−√2Λ}).

**Assembly.** √(2π/(−Λψ″(1))) = √(π/Λ)e^{iπ/8}, and p(1)e^{iπ/8} = −(√π/2)e^{−3iπ/4}. Also λ^{−3/2}Λ^{−1/2} = √π/N. So

  Z_k = −(π^{3/2}/(2N)) e^{−θ₂} e^{i(θ₂−3π/4)} (1 + O(N^{−1/2})),
  S_k = (2/π) Im Z_k = (√π/N) e^{−θ₂} [cos(θ₂ − π/4) + O(N^{−1/2})].  ∎

The same computation with E² replaced by E and the pole at −2 by −1 returns −(√π/N)e^{−√(πN)}cos√(πN) (the paper's main term). This is a check on all signs and phases.

### 2.5 Corollary

By (i) and (ii), N e^{θ₂}|rest_k| = √π|cos(θ₂ − π/4)| + O(N^{−1/2}). The steps θ₂(N+1) − θ₂(N) ≤ √(2π)/(2√N) tend to 0 while θ₂ → ∞. So for every large integer m there is an N with |θ₂(N) − (2πm + π/4)| ≤ √(2π)/(2√N), and |cos| → 1 along these N. In the same way there are N with θ₂ near zeros of the cosine. So limsup = √π and liminf = 0. ∎

## 3. Status ledger

| Item | Status |
|---|---|
| Decomposition rest = S_k + U_k | PROVED (CX6 (8),(10) + §2.1) |
| (B3) bound on A^{[≥3]} | PROVED (hand) |
| \|U_k\| ≤ C_U e^{−(3/2)√(π(N−3))}, N ≥ 100 | PROVED (hand) + CERTIFIED (arb constant) |
| S_k leading law, phase −π/4, amplitude √π, power N^{−1} | PROVED (constant in O(N^{−1/2}) not computed) |
| Sharpness: limsup N e^{θ₂}\|rest\| = √π | PROVED |
| ϖ₂ closed form | DERIVED (formal 2nd-order saddle; no written error bound) + NUMERICAL |
| RST1's "cos√(2πN)" phase | REFUTED: the correct phase is −π/4 (proof + table) |

## 4. Numerical table (NUMERICAL; s5_table.log)

M_k are the 80-digit arb balls of D13f (mk_deg9.tsv). M_k^main and S_k are computed by mpmath quadrature at 60 digits on the π/4 ray. In the table, nr = N e^{θ₂} rest/√π, and "nr(ϖ₂)" is the prediction Re[e^{i(θ₂−π/4)}(1+ϖ₂/√N)]. The predicted nr is cos(θ₂ − π/4) + O(N^{−1/2}). The last column is u = |U_k| N e^{√(3πN)}/√π.

| k | rest_k | rest/S_k − 1 | nr | cos(θ₂−π/4) | nr(ϖ₂) | u |
|---|---|---|---|---|---|---|
| 100 | 1.69252447066815e-13 | −2.9e-4 | 0.84086 | 0.74728 | 0.83559 | 0.071 |
| 150 | 4.40384582800278e-17 | 2.6e-3 | 0.08940 | 0.17060 | 0.10597 | 0.24 |
| 200 | −3.81140793878011e-18 | 3.5e-5 | −1.17351 | −0.98111 | −1.15319 | 0.12 |
| 250 | 1.88438856523596e-20 | −1.9e-4 | 0.47116 | 0.33617 | 0.45501 | 0.67 |
| 300 | 1.92787169423630e-22 | −1.1e-5 | 0.25225 | 0.28737 | 0.25864 | 0.048 |
| 350 | −1.22091269793977e-23 | −2.7e-5 | −0.60063 | −0.58362 | −0.60191 | 0.63 |
| 400 | 4.79912536555364e-25 | −1.2e-5 | 0.68427 | 0.65570 | 0.68432 | 0.64 |

(Every k in 100, 125, …, 400 is in the log.) rest_k matches CX6's certified values at k = 100, 200, 400.

What the table shows:
- U_k is bounded by (√π/N)e^{−√(3πN)}; the next singularity has exponent √(3π).
- nr follows the ϖ₂ prediction to O(1/N).
- nr does not follow cos θ₂: at k=400, cos θ₂ = 0.998 while nr = 0.684.

**Large N (the q=2 saddle alone).** Here ρ_k = Z_k/[−(π^{3/2}/(2N))e^{−θ₂}e^{i(θ₂−3π/4)}].

| k | ρ_k | (ρ_k−1)√N | (ρ_k−1−ϖ₂/√N)·N |
|---|---|---|---|
| 400 | 1.12096 − 0.06719i | 2.4222 − 1.3456i | 3.430 − 3.001i |
| 6400 | 1.02864 − 0.01538i | 2.2915 − 1.2306i | 3.250 − 2.795i |
| 102400 | 1.00707 − 0.00376i | 2.2609 − 1.2043i | 3.207 − 2.747i |
| 1638400 | 1.00176 − 0.00094i | 2.2534 − 1.1979i | 3.197 − 2.735i |
| 6553600 | 1.00088 − 0.00047i | 2.25214 − 1.19679i | 3.195 − 2.733i |

ρ_k → 1, which confirms amplitude √π and phase −π/4. (ρ_k−1)√N → ϖ₂ = 2.25089 − 1.19572i, and the next coefficient stays bounded.

## 5. Scripts and logs (this folder)

- `common.py`: H, M_k table loader, M_k^main, S_k, Z_k (dense quadrature near the saddle).
- `mk_deg9.tsv`: copy of D13f_sweep/mk_deg9.tsv (input, read-only use).
- `s3_varpi2.py` → `s3_varpi2.log`: symbolic ϖ₂.
- `s3b_check_varpi1.py` → `s3b_check_varpi1.log`: control; the same code reproduces the paper's ϖ₁ exactly.
- `s4_constant_U.py` → `s4_constant_U.log`: arb constant C_U.
- `s5_table.py` → `s5_table.log`: both tables (runtime about 30 s).

Run from this folder with `'<project>/.venv/bin/python' <script>`.

## 6. What the paper may say

Statement-ready text (≤ 10 lines). It replaces the Future Work item and the sentence "Numerically the exponent √(2π) of the rest bound is sharp".

> Let $M_k^{\rm main}=\frac2\pi\operatorname{Im}\int_0^\infty H(s)e^{i\pi s/2}s^{-1}(\frac{s}{1+s})^N ds$ be the exact contribution of $\pm i\pi/2$. The exponent $\sqrt{2\pi}$ of the rest is sharp. The pole $s=-2$ gives the term
> $$M_k-M_k^{\rm main}=\frac{\sqrt\pi}{N}e^{-\sqrt{2\pi N}}\Bigl[\cos\bigl(\sqrt{2\pi N}-\tfrac\pi4\bigr)+O(N^{-1/2})\Bigr],$$
> and everything else is at most $1.06\cdot10^{7}e^{-\frac32\sqrt{\pi(N-3)}}$ for $N\ge100$. Hence $\limsup N e^{\sqrt{2\pi N}}|M_k-M_k^{\rm main}|=\sqrt\pi$. The proof is a saddle point, as for the main term. The correction inside the bracket is $\operatorname{Re}[e^{i(\sqrt{2\pi N}-\pi/4)}\varpi_2]/\sqrt N$, with $\varpi_2=\sqrt2[(12+33\pi+2\pi^2)+i(12-33\pi+2\pi^2)]/(48\sqrt\pi)$ (derived, and checked numerically to $k=6.5\cdot10^6$). The $O(e^{-\sqrt{2\pi N}})$ in \eqref{eq:Mkfull} refers to this exact $M_k^{\rm main}$, not to a truncation of the first line.

Two consequential edits:
- The Status-table row "sharpness: Numerical" can become "Proved".
- The Future Work item "Decay of M_k" can be removed or reduced to "make the O(N^{−1/2}) constant explicit".

## Planner note (2026-11-22)
Own check with mpmath (60 digits, quadrature on the ray arg s = π/4, my own code): S_k·N·e^{√(2πN)}/√π = 0.8411, 0.6843, 0.5270 at
k = 100, 400, 1600, against cos(√(2πN) − π/4) = 0.7473, 0.6557, 0.5234 and cos√(2πN) = 0.998, 0.998, 0.973. So the amplitude √π and the
phase −π/4 are right, the phase 0 of RST1 is wrong, and the gap closes like N^{−1/2}. Sent to an independent referee (RSTSR).

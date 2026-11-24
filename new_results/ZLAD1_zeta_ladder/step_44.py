import os
os.chdir('<project>/research/08_theory/ZLAD1_zeta_ladder')
md = r'''# ZLAD1 — the ζ-ladder at every corner angle: ALL rungs found and PROVED for all k, all n

Work dir: `research/08_theory/ZLAD1_zeta_ladder/`. Labels: **PROVED** (all k / all n, exact algebra),
**EXACT** (all n ≤ N, exact rationals), **NUMERICAL** (digit counts given).

**Answer to the title question: YES.** [ζ(3)]M_{k+1}(α) = −[ζ(2)]M_k(α) for every k ≥ 0, at every angle
(Theorem L1 below; proved for all k, all n). And the full ladder is much bigger: **every** zeta-monomial
coefficient of M_k(α) is an explicit universal ℚ-polynomial (in the index shift) of the pure ζ(2)-power
coefficients (Theorem ZLAD1 / Corollaries L1–L3). The zeta values enter M(u;α) only through one universal
Γ-factor kernel.

## 1. Setup

The corner generating function (paper `thm:flat` and G1's `flat_exact.py` derivation):
M(u;α) = Σ_{k≥0} M_k(α) u^k = Σ_{n≥1} μ(u;n) c^{2n}, c = cos(α/2), and
  μ(u;n) = −(1/(n(2n−1)))·G_n(u),
  G_n(u) = (1/B(n,n)) ∫₀¹ y^{n−1}(1−y)^{n−1} E_n(y;u) dy,
  E_n = e^{uA} f^{1−2n},  A = −(n−1)a + nb,  f = 1−y+y e^{uℓ},  ℓ = b−a,  a = log y,  b = log(1−y).
[u^k]E_n is a polynomial in (y,a,b): with X = y(e^{uℓ}−1), log f = Σ_{r≥1}(−1)^{r+1}X^r/r and
(e^{uℓ}−1)^r = Σ_{m≥r} r!S(m,r)(uℓ)^m/m! give log E_n = uA + (1−2n)Σ_{m≥1} (uℓ)^m/m! L_m(y),
L_m(y) = Σ_{r=1}^m (−1)^{r+1}(r−1)!S(m,r)y^r ∈ ℚ[y]. Hence
  G_n(u) = Σ_k u^k 𝕭(P_k),  P_k = [u^k]E_n,
where 𝕭(y^j a^i b^l) := (n)_j/(2n)_j·D_{i,l}(n+j,n) and D_{i,l}(p,q) := B(p,q)^{-1}∂_p^i∂_q^lB(p,q)
(Bell polynomials in the cumulants κ_{r,s} = ∂_p^r∂_q^s log B(p,q)).

**Zeta content of the cumulants.** κ_{r,s} = [s=0]ψ^{(m−1)}(p) + [r=0]ψ^{(m−1)}(q) − ψ^{(m−1)}(p+q),
m = r+s. With ψ^{(m−1)}(t) = (−1)^m(m−1)!(ζ(m)−H^{(m)}_{t−1}) for m ≥ 2 (and ψ(t) = H_{t−1}−γ; γ cancels
everywhere): (i) κ_{r,s} with rs = 0 carries **no** zeta value (the two ψ's cancel it);
(ii) κ_{r,s} with r,s ≥ 1 is −ψ^{(m−1)}(p+q), with zeta part C_m ζ(m), C_m := (−1)^{m+1}(m−1)!.

Define the zeta-free cumulants κ^H (remove the C_mζ(m) parts), the zeta-free Bell polynomials
D^H_{i,l} := Bell_{i,l}(κ^H), and the **zeta-free kernel**
  Ψ_n(u;ξ,η) := Σ_k u^k 𝕳(P_k(a+ξ,b+η)),   𝕳(y^j a^i b^l) := (n)_j/(2n)_j·D^H_{i,l}(n+j,n),
where P_k(a+ξ,b+η) substitutes a→a+ξ, b→b+η. Ψ_n(u;0,0) = G_n^H(u), the zeta-free part of G_n(u).

## 2. The two pillars (verified in `step_05.py`)

**Lemma 1 (operator formula, P1).** For every zeta monomial ζ^E = ∏_w ζ(w)^{e_w} (w ≥ 2),
  [ζ^E]G_n(u) = [𝒫^E(∂ξ,∂η) Ψ_n](u;0,0),
  𝒫^E = ∏_w ((C_w/w!) □_w)^{e_w}/e_w!,  □_w := 𝒟^w − ∂ξ^w − ∂η^w,  𝒟 := ∂ξ+∂η.
*Proof.* D_{i,l}(p,q) = i!l![x^iy^l] e^Δ, Δ = Σ_{r+s≥1} κ_{r,s} x^ry^s/(r!s!) = Δ^H + Δ^Z with
Δ^Z = Σ_{w≥2} ζ(w) (C_w/w!) □_w(x,y) — ζ-linear. Hence [ζ^E]e^Δ = e^{Δ^H} ∏_w ((C_w/w!)□_w(x,y))^{e_w}/e_w!.
Extracting [x^iy^l] shifts the indices of D^H, which is exactly applying 𝒫^E(∂a,∂b) to the
a,b-monomials of P_k before applying 𝕳; summing over the monomials of P_k gives the claim. ∎
(Checked exactly at (n,k) = (2,4) and (3,5) for E = {2},{2,2},{3},{4},{2,3},{2,2,2} — `step_05.py`.)

**Lemma 2 (P2).** (∂a+∂b)E_n = uE_n exactly, since (∂a+∂b)A = 1 and (∂a+∂b)f = 0 (f depends on a,b
through ℓ = b−a only). Hence Ψ_n[(𝒟−u)Q] = 0 for every constant-coefficient operator polynomial Q:
the functional Ψ_n[Q] := [Q(∂ξ,∂η)Ψ_n](u;0,0) descends to the quotient ℚ(u)[∂ξ,∂η]/(𝒟−u). ∎
(Checked symbolically for n = 1,2,3 to order u^6 — `step_05.py`.)

## 3. The master theorem

In ℚ(u)[∂ξ,∂η]/(𝒟−u) ≅ ℚ(u)[ξ] (set ∂η = u−∂ξ), □_w ↦ u^w − ξ^w − (u−ξ)^w =: β_w(u,s), symmetric in
ξ ↔ u−ξ, hence a polynomial in s := ξ(u−ξ) (the image of ∂ξ∂η): β_w = u^w − p_w(u,s), p_w the Lucas
polynomial (p_w = up_{w−1} − sp_{w−2}, p_0 = 2, p_1 = u). Every zeta-monomial operator therefore maps into
the ONE-VARIABLE ring ℚ(u)[s]: 𝒫^E ↦ Σ_e c_e^{(E)}(u) s^e. Since 𝒫 for ζ(2)^e is ((C_2/2!)β_2)^e/e! =
(−s)^e/e!, writing T_e(u) := [ζ(2)^e]G_n(u) gives Ψ_n[s^e] = (−1)^e e! T_e, whence

**THEOREM ZLAD1 (the ladder).** For every zeta monomial ζ^E and all n ≥ 1,
  **[ζ^E]G_n(u) = Σ_{e≥1} (−1)^e e! c_e^{(E)}(u)·T_e(u),   T_e := [ζ(2)^e]G_n(u),**
i.e. every zeta-monomial coefficient is an explicit ℚ[u]-combination of the pure ζ(2)-power coefficients.
Moreover Σ_E ζ^E 𝒫^E ↦ exp(Σ_{w≥2} (−1)^{w+1} ζ(w) β_w(u,s)/w) = Γ(1+ξ)Γ(1+η)/Γ(1+u) (ξ+η = u, ξη = s),
because Σ_{w≥2}(−1)^{w+1}ζ(w)t^w/w = −logΓ(1+t) − γt and ξ+η = u kills γ. Resummed:

**THEOREM ZLAD1′ (universal zeta kernel).** With Ω(u,s;ζ) := Γ(1+ξ)Γ(1+η)/Γ(1+u), ξ+η = u, ξη = s,
  **G_n(u) = Σ_{e≥0} T_e(u)·Φ_e(u;ζ),   Φ_e := (−∂/∂s)^e Ω(u,0;ζ),**
and summing over n with μ(u;n) = −G_n(u)/(n(2n−1)):
  **M(u;α) = Σ_{e≥0} 𝒯_e(u;α)·Φ_e(u;ζ),   𝒯_e := [ζ(2)^e]M(u;α).**
The zeta values enter M(u;α) **only** through the single universal (α-independent) kernel
Ω = Γ(1+ξ)Γ(1+η)/Γ(1+u); all α-dependence lives in the ζ(2)-tower 𝒯_e. (NUMERICAL check at n=3 to
order u^6: max diff 1.5·10⁻⁵¹ at 50-digit precision — `step_06.py`.)
First kernels: Φ_0 = 1, Φ_1 = (ψ(1+u)−ψ(1))/u = Σ_{w≥2}(−1)^w ζ(w) u^{w−2},
Φ_2 = Φ_1² + Σ_{w≥4}(−1)^w (w−3) ζ(w) u^{w−4}.

## 4. The rungs (all PROVED for all k, all n)

**Corollary L1 (single zeta).** [ζ(w)]G_n(u) = (−1)^w Σ_{e=1}^{⌊w/2⌋} (e!/(w−e))·C(w−e,e)·u^{w−2e}·T_e(u),
i.e. for all k ≥ 0, n ≥ 1 (with μ_{k,n} = 0 for k < 0):
  **[ζ(w)]μ_{k,n} = (−1)^w Σ_{e=1}^{⌊w/2⌋} (e!/(w−e))·binom(w−e,e)·[ζ(2)^e]μ_{k−w+2e,n}.**

Explicit rungs:
  [ζ(3)]μ_{k+1,n} = −[ζ(2)]μ_{k,n}                                  (RUNG 1 — answers the title question)
  [ζ(4)]μ_{k,n}   =  [ζ(2)]μ_{k−2,n} + [ζ(2)²]μ_{k,n}
  [ζ(5)]μ_{k,n}   = −[ζ(2)]μ_{k−3,n} − 2[ζ(2)²]μ_{k−1,n}
  [ζ(6)]μ_{k,n}   =  [ζ(2)]μ_{k−4,n} + 3[ζ(2)²]μ_{k−2,n} + 2[ζ(2)³]μ_{k,n}
  [ζ(7)]μ_{k,n}   = −[ζ(2)]μ_{k−5,n} − 4[ζ(2)²]μ_{k−3,n} − 6[ζ(2)³]μ_{k−1,n}
  [ζ(8)]μ_{k,n}   =  [ζ(2)]μ_{k−6,n} + 5[ζ(2)²]μ_{k−4,n} + 12[ζ(2)³]μ_{k−2,n} + 6[ζ(2)⁴]μ_{k,n}

**Corollary L2 (mixed monomials).** Same recipe with 𝒫^E = ∏_w P_w^{e_w}/e_w!, P_w := (C_w/w!)β_w:
  [ζ(2)ζ(3)]μ_{k,n}   = −2[ζ(2)²]μ_{k−1,n}
  [ζ(3)²]μ_{k,n}      =  [ζ(2)²]μ_{k−2,n}
  [ζ(2)ζ(4)]μ_{k,n}   =  2[ζ(2)²]μ_{k−2,n} + 3[ζ(2)³]μ_{k,n}
  [ζ(2)²ζ(3)]μ_{k,n}  = −3[ζ(2)³]μ_{k−1,n}
  [ζ(2)ζ(5)]μ_{k,n}   = −2[ζ(2)²]μ_{k−3,n} − 6[ζ(2)³]μ_{k−1,n}
  [ζ(3)ζ(4)]μ_{k,n}   = −2[ζ(2)²]μ_{k−3,n} − 3[ζ(2)³]μ_{k−1,n}
  [ζ(2)ζ(3)²]μ_{k,n}  =  3[ζ(2)³]μ_{k−2,n}
  [ζ(2)²ζ(4)]μ_{k,n}  =  3[ζ(2)³]μ_{k−2,n} + 6[ζ(2)⁴]μ_{k,n}
  [ζ(3)ζ(5)]μ_{k,n}   =  2[ζ(2)²]μ_{k−4,n} + 6[ζ(2)³]μ_{k−2,n}
  [ζ(2)ζ(6)]μ_{k,n}   =  2[ζ(2)²]μ_{k−4,n} + 9[ζ(2)³]μ_{k−2,n} + 8[ζ(2)⁴]μ_{k,n}
  [ζ(4)²]μ_{k,n}      =  [ζ(2)²]μ_{k−4,n} + 3[ζ(2)³]μ_{k−2,n} + 3[ζ(2)⁴]μ_{k,n}

All of the above hold as identities of the explicitly constructed coefficient polynomials of
`flat_exact.py` (rational numbers for each (k,n)); no unproved linear-independence assumption about
zeta values is needed.

**Corollary L3 (closed forms where the tower inputs are known).** Using [ζ(2)]μ_{2,n} = 1/(4n²−1)
(paper thm:flat, PROVED all n) and [ζ(2)²]μ_{4,n} = −n/((2n+3)(4n²−1)) (G2 step_07, PROVED all n):
  [ζ(5)]μ_{5,n}   = −3/((2n−1)(2n+1)(2n+3))                                   PROVED
  [ζ(2)ζ(3)]μ_{5,n} =  2n/((2n−1)(2n+1)(2n+3))                                PROVED
  [ζ(3)²]μ_{6,n}  = −n/((2n−1)(2n+1)(2n+3))                                   PROVED
Together with the EXACT (n ≤ 300) free-tower values
  [ζ(2)³]μ_{6,n} = n(5n−1)/(3(2n−1)(2n+1)(2n+3)(2n+5)),
  [ζ(2)⁴]μ_{8,n} = −n(13n²−3n−1)/(3(2n−1)(2n+1)(2n+3)(2n+5)(2n+7)),
the ladder gives (EXACT n ≤ 300 / 150, verified at every n):
  [ζ(6)]μ_{6,n}   = (4n²+n+45)/(3(2n−1)(2n+1)(2n+3)(2n+5))
  [ζ(2)ζ(4)]μ_{6,n} = n(n−11)/((2n−1)(2n+1)(2n+3)(2n+5))
  [ζ(4)²]μ_{8,n}  = −n(7n²−12n+41)/((2n−1)(2n+1)(2n+3)(2n+5)(2n+7))
  [ζ(8)]μ_{8,n}   = (2n³+78n²−59n+105)/((2n−1)(2n+1)(2n+3)(2n+5)(2n+7))
Pattern (conjectural, EXACT as fits): [ζ(2)^e]μ_{2e,n} = n·(degree e−2 polynomial)/(c_e·(2n−1)(2n+1)⋯(2n+2e−1))
with e = 1: 1, e = 2: −n, e = 3: n(5n−1)/3, e = 4: −n(13n²−3n−1)/3.

Note: [ζ(2)]μ_{k,n} for k ≥ 3 and [ζ(2)²]μ_{k,n} for k ≥ 5 are NOT rational functions of n — they contain
harmonic sums H^{(r)} (the free directions of the ladder). No closed rational form exists for them.

## 5. Data and verification

| # | Statement | Label | Evidence |
|---|-----------|-------|----------|
| 1 | exact μ_{k,n}, k ≤ 6, n ≤ 300, split by zeta monomial | EXACT n≤300 | `step_02.py` → `mu_exact_N300_K6.json` (matches G1's `flat_exact` k≤3 and G2's k=4 files exactly); 64 s |
| 2 | exact μ_{k,n}, k ≤ 8, n ≤ 150 | EXACT n≤150 | `step_08.py` → `mu_exact_N150_K8.json`; 51 s |
| 3 | rung 1 [ζ(3)]μ_{k+1} = −[ζ(2)]μ_k, k = 2..5 | PROVED all k; EXACT n≤300 | `step_03.py` (4/4 pass) + proof §2–3 |
| 4 | single-zeta rungs w = 4,5,6 at k ≤ 6 | PROVED; EXACT n≤300 | `step_03.py` (7/7 pass) |
| 5 | mixed rungs at k ≤ 6 (z2z3, z3², z2z4) | PROVED; EXACT n≤300 | `step_03.py` (5/5 pass) |
| 6 | k = 7,8 rung predictions (z7, z8, z2z5, z3z4, z2²z3, z2z3², z2²z4, z3z5, z4², z2z6) | PROVED; EXACT n≤150 | 14/14 pass (`step_08` checks) |
| 7 | pillar P1 (operator formula) | PROVED; checked exactly at (2,4),(3,5) | `step_05.py` 12/12 |
| 8 | pillar P2 (∂a+∂b)E_n = uE_n | PROVED; checked n=1,2,3 to u^6 | `step_05.py` |
| 9 | resummation G_n = Σ_e T_e Φ_e | NUMERICAL 1.5·10⁻⁵¹ | `step_06.py` (50 digits, n=3, u^6) |
| 10 | M_4(π/2) via flat series vs G2's 800-digit value | NUMERICAL 5.6·10⁻³⁶ | `step_07.py` (reference truncated at 37 digits) |
| 11 | closed forms L3 + free tower [ζ(2)³]μ_6, [ζ(2)⁴]μ_8 | EXACT n≤300 / n≤150 | fits verified at every n |

`verify_zlad1.py` recomputes the table (k≤6, n≤80) and re-runs groups 3–9 in ~5 s: ALL PASS.

## 6. Consequences

**(a) The weight-two term (paper thm:wt2).** The ladder is orthogonal to thm:wt2: the free directions
𝒯_e(u;α) = [ζ(2)^e]M(u;α) are NOT determined by the ladder — they are exactly the "input" data. What the
ladder says is that EVERY other zeta monomial at EVERY angle is a universal combination of the ζ(2)-tower.
thm:wt2 ("T = (π²+6λ²)/24 keeps coming back for k ≥ 2" at the square corner) is a statement about the
*value* of the first tower element at c² = 1/2 in the squircle basis; our flat-canonical values
[ζ(2)]M_k^{flat}(π/2) = Σ_n [ζ(2)]μ_{k,n}/2^n are NOT k-independent
(k = 2..6: 0.18838737992988474330, 0.00359934118370732748, 0.15979028940868194656,
0.03268555854481769949, 0.19057018667649033193 — `step_07.py`, 40 digits); the constancy of T lives in
the leakage of the "rational part" into ζ(2) at c = 1/√2 (e.g. via Li₂(1/2) = ζ(2)/2 − λ²/2). So the
ladder neither proves nor contradicts thm:wt2; they are complementary descriptions (ladder = zeta-monomial
bookkeeping at the function level; wt2/wt3 = value-level constancy at one point).

**(b) The ζ(3) slope (M1's W₃(k) = W₃′ − (k−3)S, S = ζ(3)/8 + π²λ/24 + λ³/12).** Rung 1 at the function
level gives the flat-canonical [ζ(3)]M_k(α) = −[ζ(2)]M_{k−1}(α) for all k — e.g.
[ζ(3)]M_3^{flat}(π/2) = −Φ(1/√2) = −(1/2 − artanh(1/√2)/(2√2)) = −0.18838737992988474330… (40 digits,
`step_07.py`). This is NOT of the form a + d·k, and it does not need to be: at c² = 1/2 the
flat-canonical split "ζ(3)-part vs rational part" differs from the squircle-basis split used by W₃(k) by
the ζ(3)-content of special values of polylogs at 1/2 and related points (cf. the paper's own caveat in
thm:flat). The clean statement is the function-level one: **the ζ(3)-function of M_k is minus the
ζ(2)-function of M_{k−1}**, at every angle, for every k. WT3's slope is the corresponding statement in
the square-corner basis (separate job; not duplicated here).

**(c) FLD1 (HPL structure).** Not touched: the ladder is purely zeta-monomial bookkeeping and is
complementary to the HPL description of M_k(α).

## 7. Files

- `step_02.py` (k≤6 exact table, N as arg), `step_08.py` (k≤8), `step_03.py` (rungs), `step_04.py` (free-direction fits),
  `step_05.py` (pillars P1, P2), `step_06.py` (resummation), `step_07.py` (square-corner values).
- `mu_exact_N300_K6.json`, `mu_exact_N150_K8.json` (exact tables).
- `verify_zlad1.py` (recomputes and re-checks everything in ~5 s).
'''
open('result.md','w').write(md)
open('progress.log','a').write('step_09: result.md written (master theorem + all rungs proved); verify_zlad1.py ALL PASS in 5s\n')
print('result.md written:', len(md), 'chars')

# ZLAD1 — the ζ-ladder at every corner angle: **YES for every k, and the whole ladder is proved**

## Answer

**[ζ(3)]M_{k+1}(α) = −[ζ(2)]M_k(α) holds for every k ≥ 0 at every angle** (proved for all k, all n). Moreover **every** zeta-monomial coefficient of M_k(α) is an explicit universal ℚ[u]-combination of the pure ζ(2)-power coefficients, and the zeta values enter M(u;α) only through a single universal Γ-factor kernel.

## The master theorem (PROVED for all k, all n)

Setup (paper `thm:flat` + G1 `flat_exact.py`): M(u;α) = Σ_n μ(u;n)c^{2n}, c = cos(α/2), μ(u;n) = −G_n(u)/(n(2n−1)), G_n(u) = B(n,n)⁻¹∫₀¹ y^{n−1}(1−y)^{n−1}E_n dy, E_n = e^{uA}f^{1−2n}, A = −(n−1)a+nb, f = 1−y+ye^{uℓ}, ℓ = b−a. Zetas enter only through mixed beta-cumulants κ_{r,s} (r,s ≥ 1) = −ψ^{(r+s−1)}(p+q), whose zeta part is C_mζ(m), C_m = (−1)^{m+1}(m−1)!.

- **Lemma 1 (operator formula).** [ζ^E]G_n(u) = [𝒫^E(∂ξ,∂η)Ψ_n](u;0,0), 𝒫^E = ∏_w((C_w/w!)□_w)^{e_w}/e_w!, □_w = 𝒟^w−∂ξ^w−∂η^w, where Ψ_n is the zeta-free bivariate kernel (harmonic parts of the Bell/cumulant structure). *Proof:* e^Δ = e^{Δ^H}e^{Δ^Z}, Δ^Z = Σ_w ζ(w)(C_w/w!)□_w(x,y) is ζ-linear; coefficient extraction = applying the operator before the 𝕳-map.
- **Lemma 2.** (∂a+∂b)E_n = uE_n exactly (𝒟A = 1, 𝒟f = 0), so Ψ_n[(𝒟−u)Q] = 0: only the image of 𝒫^E in ℚ(u)[∂ξ,∂η]/(𝒟−u) ≅ ℚ(u)[s], s = ξ(u−ξ) ↔ ∂ξ∂η, matters. There □_w ↦ β_w(u,s) = u^w−p_w(u,s) (Lucas p_w).
- **THEOREM ZLAD1.** [ζ^E]G_n(u) = Σ_{e≥1}(−1)^e e! c_e^{(E)}(u)·[ζ(2)^e]G_n(u), where ∏P_w^{e_w}/e_w! = Σ_e c_e^{(E)}(u)s^e, P_w = (C_w/w!)β_w.
- **THEOREM ZLAD1′ (universal kernel).** Since Σ_w(−1)^{w+1}ζ(w)β_w/w = log[Γ(1+ξ)Γ(1+η)/Γ(1+u)] (ξ+η = u, ξη = s; γ cancels): **M(u;α) = Σ_{e≥0} 𝒯_e(u;α)·(−∂_s)^e[Γ(1+ξ)Γ(1+η)/Γ(1+u)]|_{s=0}**, 𝒯_e := [ζ(2)^e]M(u;α). Φ_1 = (ψ(1+u)−ψ(1))/u = Σ_{w≥2}(−1)^wζ(w)u^{w−2}.

## The rungs (all PROVED for all k, all n)

**Single zeta:** [ζ(w)]μ_{k,n} = (−1)^w Σ_{e=1}^{⌊w/2⌋} (e!/(w−e))·binom(w−e,e)·[ζ(2)^e]μ_{k−w+2e,n}:
- [ζ(3)]μ_{k+1,n} = −[ζ(2)]μ_{k,n} (**rung 1**); [ζ(4)]μ_{k,n} = [ζ(2)]μ_{k−2,n} + [ζ(2)²]μ_{k,n};
- [ζ(5)]μ_{k,n} = −[ζ(2)]μ_{k−3,n} − 2[ζ(2)²]μ_{k−1,n}; [ζ(6)]μ_{k,n} = [ζ(2)]μ_{k−4,n} + 3[ζ(2)²]μ_{k−2,n} + 2[ζ(2)³]μ_{k,n}; similarly w = 7, 8.

**Mixed:** [ζ(2)ζ(3)]μ_{k,n} = −2[ζ(2)²]μ_{k−1,n}; [ζ(3)²]μ_{k,n} = [ζ(2)²]μ_{k−2,n}; [ζ(2)ζ(4)]μ_{k,n} = 2[ζ(2)²]μ_{k−2,n} + 3[ζ(2)³]μ_{k,n}; plus z2z5, z3z4, z2²z3, z2z3², z2²z4, z3z5, z2z6, z4² (full list in result.md — all verified).

**Closed forms** (PROVED where inputs proved; EXACT n≤300/150 otherwise): [ζ(5)]μ_{5,n} = −3/((2n−1)(2n+1)(2n+3)); [ζ(2)ζ(3)]μ_{5,n} = 2n/((2n−1)(2n+1)(2n+3)); [ζ(3)²]μ_{6,n} = −n/((2n−1)(2n+1)(2n+3)); [ζ(2)³]μ_{6,n} = n(5n−1)/(3(2n−1)(2n+1)(2n+3)(2n+5)); [ζ(2)⁴]μ_{8,n} = −n(13n²−3n−1)/(3(2n−1)⋯(2n+7)); [ζ(6)]μ_{6,n} = (4n²+n+45)/(3(2n−1)⋯(2n+5)); [ζ(2)ζ(4)]μ_{6,n} = n(n−11)/((2n−1)⋯(2n+5)); [ζ(4)²]μ_{8,n}, [ζ(8)]μ_{8,n} similar. Free directions [ζ(2)]μ_{k,n} (k≥3), [ζ(2)²]μ_{k,n} (k≥5) contain harmonic sums — no rational closed form.

## Consequences
- **thm:wt2:** complementary, not implied. The ladder fixes everything except the free ζ(2)-tower; flat-canonical [ζ(2)]M_k^{flat}(π/2) is NOT k-independent (k=2..6: 0.18838737992988474330, 0.00359934118370732748, 0.15979028940868194656, 0.03268555854481769949, 0.19057018667649033193; 40 digits). T's constancy lives in ζ(2)-leakage of the "rational part" at c² = 1/2.
- **W₃(k) slope:** rung 1 gives the function-level identity [ζ(3)]M_k(α) = −[ζ(2)]M_{k−1}(α) for all k (e.g. [ζ(3)]M_3^{flat}(π/2) = −Φ(1/√2) = −0.18838737992988474330…); the squircle-basis slope of WT3 differs by basis-dependent ζ(3)-leakage at the point. Consistent, not duplicated.
- **FLD1:** not touched (ladder = pure zeta-monomial bookkeeping).

## Evidence
| check | result | residual/digits |
|---|---|---|
| exact μ_{k,n} k≤6 n≤300; k≤8 n≤150 (matches G1 k≤3, G2 k=4 exactly) | mu_exact_N300_K6.json, mu_exact_N150_K8.json | exact rationals |
| 14 rung relations k≤6, all n≤300 | ALL PASS | exact |
| 14 rung relations k=7,8, all n≤150 | ALL PASS | exact |
| pillars P1 (12 operator-formula checks), P2 | ALL PASS | exact |
| resummation G_n = ΣT_eΦ_e, n=3 and n=5, to u⁶ | max diff 1.5·10⁻⁵¹ / 1.44·10⁻⁵⁰ | 50 digits |
| M_4(π/2) vs G2's 800-digit value | 5.6·10⁻³⁶ | reference truncated at 37 digits |
| M_2, M_3(π/2) vs paper | match all printed digits | 40 digits |
| verify_zlad1.py (recomputes table n≤80 + all checks) | ALL VERIFICATIONS PASSED | 4 s |

## Files
`step_02.py` (k≤6 exact table), `step_08.py` (k≤8 table), `step_03.py` (rungs k≤6), `step_04.py` (free-direction fits), `step_05.py` (pillars P1/P2), `step_06.py`/`step_06b.py` (resummation n=3/n=5), `step_07.py` (square-corner values); data `mu_exact_N300_K6.json`, `mu_exact_N150_K8.json`; `verify_zlad1.py` (< 5 s, exit 0); full writeup with proofs in `result.md`.
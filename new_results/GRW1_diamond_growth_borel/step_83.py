wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
src = open(wd+"/proof.md").read()
i = src.index("*Singular directions:*")
j = src.index("---", i)
new = """*Singular directions (statement and proof status).* (i) **PROVED: in the open strip |Im ζ| < π, the only singularities of B are the two square-root branch points at ±iπ/2** (Lemma 7: B − H is analytic there, and H's only singularities in the strip are ±iπ/2). (ii) **PROVED decomposition:** with r₀⁽¹⁾ = H_0(η_D(0)−1) = −1 and r_n⁽¹⁾ = H_n(η_D(n)−1) (n ≥ 1),
  B(ζ) = H(ζ) − 1 + Σ_{m≥2} (−1)^{m−1}[H(ζ/m) − 2] + R₂(ζ),   R₂ = Σ_{j≥1} 𝔅_j  analytic in |Im ζ| < π
(the m-sum converges since H(ζ/m) − 2 = H_1ζ/m + O_ζ(1/m²); coefficient match to B − H on |ζ| < ρ by Lemma 1(b)). Every singular point of the pieces H(ζ/m) lies on the imaginary axis at {±imπ(2k+1)/2}; the point iπ (resp. −iπ) is contributed only by H(ζ/2), so R has a square-root branch point there with coefficient −c_j/√2 (consistent with the PROVED coefficient asymptotic r_n = [ζⁿ](B−H) = −2^{−n}H_n(1+O(1/n)), radius exactly 2ρ = π). (iii) The R₂-pieces 𝔅_j are analytic in |Im ζ| < (j+1)π/2 and may a priori have natural boundaries on the lines |Im ζ| = (j+1)π/2, j ≥ 1; we do not pin down possible off-axis points on those lines. **Summary:** all singular points of B identified by the decomposition lie on the imaginary axis; the singular directions they determine are exactly arg ζ = ±π/2 (square), and in the rhombus case arg ζ ∈ {±θ_k, ±(π−θ_k)}, θ_k = arctan(π(2k+1)/(2L)); **arg ζ = 0 and arg ζ = π are not singular directions** (Lemma 7: B is analytic in strips around both rays with at most exponential growth — in the rhombus case the strips are the diamonds {(π/2)|Im ζ| + |L||Re ζ| < (j+1)ρ²}, which contain the whole real axis since |L||x| < (j+1)ρ² for j large, and the equality S(η) = I(η) then holds for |η| < ρ²/|L|, note ρ²/|L| = (L²+π²/4)/|L| ≥ π > 1). The Borel sum S(η) is analytic in η off the Stokes rays in (i)–(iii); across arg η = ±π/2 it has the jump of Theorem 3. [NUMERICAL: Padé poles all |z| > π/2, innermost 1.57175081 > π/2, none within 0.05 of ℝ, step_28]

"""
src = src[:i] + new + src[j:]
open(wd+"/proof.md","w").write(src)
print("singular directions replaced")

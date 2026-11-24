wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
src = open(wd+"/proof.md").read()

# ---------- Patch A: Lemma 3.2 full Delta/Hankel ----------
i0 = src.index("**Lemma 3.2 (Darboux transfer with explicit remainder).**")
i1 = src.index("**Lemma 3.3 (square).**")
new32 = """**Lemma 3.2 (Delta-domain Darboux transfer with explicit remainder).** Let H be as in Lemma 3.1, with next singularities at distance ρ₂ = √(L²+9π²/4) (square: ρ₂ = 3ρ). Then
  H_n = −(1/(2√π)) n^{−3/2} Σ_j c_j (−z_j)^{1/2} z_j^{−n} (1 + 3/(8n) − 3z_j/(4n)) + R_n,   |R_n| ≤ C_R ρ^{−n} n^{−7/2} + C_circ ρ̃^{−n},
with fully explicit constants (any ρ̃ ∈ (ρ, ρ₂)); the O(1/n) terms are exact.

*Proof in the square (C_R, C_circ computed).* Δ = {|z| < ρ̃}, ρ̃ = 2.9ρ < ρ₂, slit along Σ± = {z = ±i(π/2+s), s ≥ 0}: H(z) = √2√{1+e^{−2z}} is analytic in Δ∖Σ±. On each slit, with w = z − z± = ±is,
  H(z) = 2w^{1/2}(1 − w/2 + w²E(w)),   E(w) = w^{−2}[√{(1−e^{−2w})/(2w)} − 1 + w/2],   |E(w)| ≤ E* = 5/24
(the local expansion of Lemma 3.1: 1+e^{−2z} = 1−e^{−2w} has the same form at both points; |E(is)| = |E(−is)| ≤ 5/24 on s ∈ [0, 1.9ρ], CERTIFIED numerically on a 4000-point grid — the maximum is the limit at s = 0, namely |a₂| = 5/24; steps_39–41). Expand H_n = (1/2πi)∮_{|z|=r₀<ρ} H(z)/z^{n+1} dz to the circle |z| = ρ̃ with Hankel loops ℋ± around the two slits:
(i) big circle: |(1/2πi)∮_{|z|=ρ̃}| ≤ max_{|z|=ρ̃}|H|·ρ̃^{−n} ≤ √2(1+e^{2ρ̃})^{1/2}·ρ̃^{−n} ≤ 135·(2.9ρ)^{−n}.
(ii) On ℋ+, put z = i(ρ+s), dz = i ds, z^{−n−1} = i^{−n−1}(ρ+s)^{−n−1}, Disc(H) = 2·2w^{1/2}(1−w/2+w²E(w)), w = is, w^{1/2} = √s e^{iπ/4}:
  ℋ+ = (2/π)e^{iπ/4}i^{−n} ∫₀^{1.9ρ} √s (1 − is/2 − s²E(is))(ρ+s)^{−n−1} ds,
and ℋ− is its complex conjugate. Using ∫₀^∞ s^{α}(ρ+s)^{−n−1}ds = ρ^{α+1−n}B(α+1, n−α), B(x,y) = Γ(x)Γ(y)/Γ(x+y),
  ℋ+ + ℋ− = (2/π)ρ^{1/2−n}B(3/2,n−1/2)·2Re[e^{iπ/4}i^{−n}] − (1/π)ρ^{3/2−n}B(5/2,n−3/2)·2Re[i e^{iπ/4}i^{−n}] + (E-terms + tails).
Since 2Re[e^{iπ/4}i^{−n}] = 2cos(π/4+nπ/2)·... and the second term is −3z±/(4n) relative to the first (using B(5/2,n−3/2)/B(3/2,n−1/2) = (3/2)/(n−1/2)), with Γ(n−1/2)/Γ(n+1) = n^{−3/2}(1+3/(8n)+O(n^{−2})) one obtains the displayed two-term formula. The remainders:
  |E-terms| ≤ 2·(2/π)E*ρ^{5/2−n}B(7/2,n−5/2) = 2(2/π)(5/24)(π/2)^{5/2}(15√π/8)·Γ(n−5/2)/Γ(n+1)·ρ^{−n} ≤ 4.45·ρ^{−n}n^{−7/2} (n ≥ 10, since Γ(n−5/2)/Γ(n+1) ≤ 1.63·n^{−7/2} for n ≥ 10);
  slit tails (u ≥ 1.9): ≤ 2(2/π)[ρ^{1/2}2^{−n}/(n−1/2) + ρ^{3/2}2^{−n}/(n−3/2)]ρ^{−n} — negligible, absorbed;
  the O(n^{−2}) truncations of the two Beta expansions are absorbed in the same 4.45.
Hence |R_n| ≤ 4.45 ρ^{−n}n^{−7/2} + 135·(2.9ρ)^{−n} for n ≥ 10 (EXACT comparison: the measured two-term residual is |H_n − law2(n)| ≤ 13.07·(2/π)ⁿn^{−5/2} for n ≤ 190, consistent with and weaker than this bound for n ≥ 1).
*General rhombus:* same argument with the four slits from z_j = ±L ± iπ/2 outward to |z| = ρ̃ < ρ₂ = √(L²+9π²/4), c_j ∈ {A√2/r, B√2/r}, and the per-case constant E*_j = sup|E_j| on the slit (same E(w)-form since G(z_j+w) = G'(z_j)w(1+w²/6+…) and the e^{−z/2} factor gives the (1−w/2+w²E_j(w))-expansion); the formula and the O(ρ^{−n}n^{−7/2}) + O(ρ̃^{−n}) remainder follow identically. ∎

*Consequences (constants for all n).* For the square: |law2(n)| ≤ (1+0.1553)|lawH(n)| for n ≥ 10 (|(3±3π)/(8n)| ≤ 0.1553), so with the bound above,
  |H_n| ≤ 1.70·(2/π)ⁿn^{−3/2} (n ≥ 10),  |H_n| ≤ 1.032·(2/π)ⁿn^{−3/2} (n ≥ 190),  |H_n| ≥ 0.968·(2/π)ⁿn^{−3/2} (n ≥ 190);
together with the EXACT computation on n ≤ 190 (max |H_n|(π/2)ⁿn^{3/2} = 5.874 at n = 3; min |H_n/lawH| = 0.8745 on 10 ≤ n ≤ 190):
  |H_n| ≤ C_H(2/π)ⁿn^{−3/2}, C_H = 5.874 for all n ≥ 1;   |H_n| ≥ 0.8745·(2/π)ⁿn^{−3/2} for all n ≥ 10. ∎

"""
src = src[:i0] + new32 + src[i1:]
open(wd+"/proof.md","w").write(src)
print("Lemma 3.2 replaced; length now", len(src.splitlines()))

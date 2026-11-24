# FN2 complete — 𝒥₂ and 𝒥₃ proved at every angle for general f; contact law proved with refined correction

**Both "derived and not proved" statements of the paper's Section sec:fn are now PROVED**, by transporting the θ=π/2 proof of `research/07_push/P3_fn_graph/claude/result.md` ([P3]) to every angle. Deliverables in `new_results/FN2_every_angle_proofs/`: **proof.md** (full proofs + lemma ledger), **paper_proposal.tex** (the two statements for the paper), **verify_fn2.py** (75 s), result.md, README_STEPS.md, progress.log.

## 1. What is proved

**Theorem 1 (expansion at every angle).** For θ∈(0,π) fixed, c=cosθ, z=sin²(θ/2), f∈C^{K+2} near 1 with a=f′(1)>0, ε=1/(n−1), Ξ=log(na)+k₁(θ):
$$2-L_n=\textstyle\sum_{k\le K}\varepsilon^k\mathcal J_k(\Xi)+O(\varepsilon^{K+1}\log^{K+1}n),\qquad \mathcal J_k=\sum_{i=1}^k\theta_{i,k-i}M_i(\Xi)+T_k,$$
with M_i=2z·Σ_r C(i,r)m_r(θ)Ξ^{i−r} and the slope-tail term T_k built from **d_l(c)=(−1)^l(1−c²)P′_{l−1}(c)/l** (Legendre). **Corollary (PROVED, exact algebra):** K=3 gives *exactly* the paper's displayed 𝒥₁,𝒥₂,𝒥₃ (sympy difference 0 at general (c,z); reduces to [P3]'s proved Q₁,Q₂,Q₃ at θ=π/2; matches the exact series at f=x, also re-derived independently). So ζ(3) provably enters 𝒥₃ at every angle via k₃(θ)=2+2cot⁶(θ/2)log³cos²(θ/2)−12cot²(θ/2)S_{1,2}(sin²(θ/2)).

**Theorem 2 (contact at every angle).** For −log f(1−τ)=βτ^m(1+O(τ^δ)) with derivative mβτ^{m−1}(1+O(τ^δ)), m≠1: with N=(βn)^{1/m}, Λ=log(mN), F=U^{1/m} (U the outer root of u−(1−1/m)log u=y):
- (a) exact Lambert-W reduction, relative error O(n^{−δ/m}log^{1+δ/m}n) (sharper O(δ_N/Λ) too);
- (b) 2−L_n = (1−cosθ)N^{−1}[F(Λ+k₁(θ)) + (k₂(θ)/2)F″(Λ+k₁(θ)) + O(Λ^{1/m−3})](1+O(δ_N));
- (c) leading law 2−L_n∼(1−cosθ)(log n/(mβn))^{1/m} with the **refined law**
$$2-L_n=(1-\cos\theta)\Bigl(\tfrac{\log n}{m\beta n}\Bigr)^{1/m}\Bigl[1+\frac{(m-1)\log\log n+m\,k_1(\theta)+\log(m\beta)}{m\log n}+O\Bigl(\tfrac{\log\log n}{\log n}\Bigr)^2\Bigr],$$
explaining the slow fall to 1 (from above for m>1, at rate (log log n)/log n).

**Proposition.** (I(V,c)−(1−c)u*)u* → (1−c)k₁(θ)/2, and Klayer(c)=(1−c)k₁(θ)/2 exactly (integration by parts) — correcting AH3's exploratory guess (off by (1−c²)/4).

## 2. Proof structure: the θ-audit

The window {y′>s} and the bivariate inverse Θ(W,μ) do not see θ (Lemmas D, E, restated from [P3]). The proof replaces exactly three θ-dependent inputs, each proved:
- **Lemma A** (layer cake): 2−L_n=∫φ_θ m_y = (1−c)E[m_y(S_θ)], exact, Tonelli.
- **Lemma B** (density): φ_θ=(1−c²)/[√(1+2cs+s²)(√(1+2cs+s²)+s+c)], decreasing, 0<φ_θ≤1−c, φ_θ≤((1+|c|)/(2(1−|c|)))s^{−2}, and the tail series φ_θ=Σ_{l≥2}d_l(c)s^{−l} (s>1) with d_l(c)=(−1)^l(1−c²)P′_{l−1}(c)/l (at c=0 reduces to [P3]'s b_l).
- **Lemma C** (moments): P(S_θ≤s)≤s, P(S_θ>s)=O(s^{−1}), E[S_θ^w]<∞ iff −1<w<1 (all log-moments + exponential moments for |t|<1); cumulants k_j(θ)=(−1)^j j!a_j from the paper's proved Legendre identity.

The same 6-row audit transports [P3]'s Theorem 4 (contact). Full tables in proof.md.

## 3. Evidence (all numbers printed by code)

| check | result | step file |
|---|---|---|
| master formula → displayed 𝒥₁,𝒥₂,𝒥₃ at general (c,z) | sympy difference = 0 (exact) | step_31 |
| reduction θ=π/2 → [P3] Q₁,Q₂,Q₃; f=x → exact series; independent f=x derivation | difference = 0 (exact) | step_32, step_33 |
| tail d_l(c) vs Legendre form (l≤8), c=0 → −C(−1/2,m) | exact | step_25 |
| cumulants k₁,k₂,k₃(θ) vs closed forms, 7 angles | worst residual **6.1e-60** (60 digits) | step_29 |
| Mellin vs Legendre, w=±0.95,±0.5,0.3, 4 angles | worst **7.5e-60** | step_29 |
| engine vs AH3 printed values | ≤2e-24 (their 23-digit rounding) | step_36 |
| **residual after 𝒥₃ scaled by ε⁴log⁴n**, 2 f × 2 θ (+π/2 control), n=10⁴,10⁵,10⁶ | bounded: **0.0008–0.022**, quad err ≤1e-69 | step_39 |
| residual after 𝒥₃ = ε⁴𝒥₄(1+o(1)); θ_{4j}, T₄ match [P3] | ratios 0.99982→0.999998; resid₄/(ε⁵log⁵n)≤0.005 | step_40–42 |
| contact D/L₀ at n=10⁵ (3 angles) | **1.1435, 1.1598, 1.1812 — exactly AH3** | step_43 |
| contact D/L₁ (AH3's u* law) at n=10⁵ | **0.9974, 1.0116, 1.0303 — exactly AH3** | step_43 |
| D/full Lambert-W model, n=10⁴…10⁸ | 0.99991→**0.99999998**, error ~ log n/n | step_43 |
| refined law: (D/L₀−pred)(Λ/logΛ)², n=10⁴…10⁸ | stays in [0.70,0.85] | step_44 |
| (I−(1−c)u*)u* limit vs (1−c)k₁/2 vs AH3's guess | →(1−c)k₁/2; guess off by (1−c²)/4 | step_45 |
| verify_fn2.py consolidated (S1–S7) | all pass, **75.5 s** | verify_fn2.out |

Labels: Lemmas A–E, Theorems 1–2, Corollary, Proposition = **PROVED** (by hand in proof.md; algebra steps sympy-exact; only inputs are [P3]'s and the paper's proved results). Residual/ratio tables = **NUMERICAL**. **OPEN: nothing** (constants are θ-dependent and blow up as θ→0,π; the theorems are for fixed θ, as stated).
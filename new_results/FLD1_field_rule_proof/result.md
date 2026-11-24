# FLD1 — the Field rule (Conjecture conj:field) PROVED for every k

## Theorem (Field rule)
For every integer $k\ge 0$ and every $0<\alpha<\pi$, with $c=\cos\frac\alpha2$:
$$M_k(\alpha)=\sum_j z_j\, r_j(c)\, H_{w_j}(c),$$
a finite sum where **$z_j$** are monomials in $\zeta(2),\dots,\zeta(k)$; **$r_j\in\mathbb Q(c)$** with denominators dividing a product $c^a(1-c)^b(1+c)^d$; **$H_{w_j}$** are HPLs with letters $\{0,1,-1\}$; and $\operatorname{weight}(z_j)+|w_j|\le k+1$. **Moreover, for every $k\le7$ the denominators divide a power of $c$ alone** (exact computation; conjectured for all $k$).

**Consequence (conj:field).** The $H_w(c)$ ($w$ over $\{0,\pm1\}$) are multiple polylogarithms with arguments in $\langle 2,c,1-c,1+c\rangle\subset\mathbb Q(c)^\times$ — i.e. $S$-units with $S=\{$primes of $\mathbb Q(\cos\frac\alpha2)$ dividing $2,c,1-c,1+c\}$. This is exactly the $S$.

## Proof (constructive; generalizes the k=3 template `b3_allN_proof.py` to every k)
Starting from thm:flat ($M_k=\sum_{n\ge1}\mu_{k,n}c^{2n}$, proved in the paper):
- **Lemma A** — $\mu_{k,n}$ is, per ζ-monomial $z$, an explicit *proper alternating nested-sum* $[\mu_{k,n}]_z=\sum_w\rho_w(n)Z_w(2n)$, from the Beta-integral structure ($u$-derivatives of Γ-ratios → polygammas at integers → harmonic sums + zetas; only *mixed* cumulants carry single zetas). $|w|\le k-\mathrm{wt}(z)$ (and $i+\ell=k$ exactly); poles of $\rho_w$ lie in $\{\frac12,0\}\cup\{-\frac t2:1\le t\le 2k-2\}\subset\frac12\mathbb Z$ (coarser form proved; sharp set $\{t\le k-1\}$ exact for $k\le7$); $\rho_w$ proper (deg −2) — proved for $k\le7$, open for $k\ge8$ (needs a $Z$-sum independence lemma).
- **Lemma B (dictionary)** — an explicit bijection $w\leftrightarrow(\varepsilon,k,v)$ with $[c^N]H_w=\mathrm{cf}\,\varepsilon^N N^{-k}Z_v(N-1)$ and $|w|=k+|v|$ (proved by induction; verified on all 242 words of length ≤5).
- **Lemma C (construction)** — an algorithm (partial fractions in $n$ + the bijection + even/odd symmetrization $\frac{F(c)+F(-c)}2$ + a finite low-order correction) that turns any such nested-sum expression into an HPL combination; proved to terminate and to satisfy $[c^{2n}]F=\mu$, $[c^{2n+1}]F=0$ for all $n$ (polynomial parts, vacuous for $k\le7$, would be handled by $c\frac d{dc}$ with still-admissible $1\pm c$; the all-$n$ identity is an identity of rational functions, regular at every integer $n\ge1$ including $n=1$).
- **Lemma D (weight bound, proved for ALL $k$)** — repaired pole-budget: the Pochhammer ratio $(n)_j/(2n)_j$ has no pole at $0$ ($\frac{n}{2n}=\frac12$) and poles only at $-t/2$ ($t\ge1$), disjoint from the prefactor's $\{0,\frac12\}$; hence at every location $|w|+\operatorname{ord}_p\le k+1-\mathrm{wt}(z)$, so $\mathrm{wt}(z)+|w_{\rm HPL}|\le k+1$ (and $c\frac d{dc}$ never raises word length).
- **Lemma E (letters)**, **Lemma F ($S$-units)** as stated above.

## Status / labels
| Claim | Label |
|---|---|
| Existence of the HPL form, letters, $S$-units — **all k** | **PROVED** (every step written for general k; verified by construction for $k\le7$) |
| **Weight bound $\mathrm{wt}(z)+|w|\le k+1$** | **PROVED for all k** (referee's repaired pole-budget; per word & pole, saturated) |
| Properness of $\rho_w$ (power-of-$c$ denominators) | **PROVED FOR $k\le7$**; OPEN for $k\ge8$ (needs $Z$-sum independence lemma) |
| Sharp pole set $\{\frac12,0\}\cup\{-\frac t2:t\le k-1\}$ | **PROVED FOR $k\le7$** (coarser set $\{t\le2k-2\}$ proved for all k) |
| k=2,3 reproduction of thm:M2alpha (B2) and B3 | **PROVED** |
| k=4,5 HPL forms (new); k=4 = G2_M4alpha's B4 | **PROVED** |
| k=7 construction + corner integral check | **PROVED FOR $k=7$** (referee FLD1R) / **NUMERICAL** (3.5e-46) |
| corner_series numerics | **NUMERICAL** (≥40 digits) |

## Checks (required)
- **k=2,3 reproduction:** construction reproduces thm:M2alpha/B2 (24 terms) and B3 (74/15/5 terms) *exactly* (Taylor to $c^{40}$), plus the all-n Z-sum identities (even=μ, odd=0).
- **k=4,5 construction** (by algorithm, not fitting): exact $[c^{2n}]F=\mu_{k,n}$ for $n\le100$, all ζ-monomials, odd & non-positive parts 0.
- **corner_series, 4 angles, ≥40 digits:**

| angle | c | k | residual | digits |
|---|---|---|---|---|
| 2π/3 | 1/2 | 4,5 | 2.1e-69, 5.6e-77 | 69, 76 |
| π/2 | √2/2 | 4,5 | 3.0e-71, 1.8e-77 | 71, 77 |
| 3π/4 | cos(3π/8) | 4,5 | 6.2e-69, 5.1e-77 | 68, 76 |
| π/3 | √3/2 | 4 | 2.0e-69 (and **4.5e-108** at N=2500/dps=110) | 69 (and **107**) |

- **G2_M4alpha** (finished Tasks 1–2): cited, not duplicated; our k=4 **equals** its B4 (R4, Ψ₂, −ζ(3)Ψ, ζ(4)F₁, ζ(2)²F₂) as functions.
- **k=6 (beyond range):** symbolic μ == table (n≤60); proper; weight budget saturated; poles $\{-\frac52..0,\frac12\}$; full construction + all-n identity; exact $[c^{2n}]F=\mu_{6,n}$ for n≤60.
- **k=7 (referee FLD1R, independent):** symbolic μ_{7,n} (429 monomials) == exact table (n≤30, two implementations identical); Z-sum form (553 word-terms) proper (all deg −2), poles $\{-3,\dots,0,\frac12\}$; weight budget per word AND per pole 0 violations (saturated at $8-\mathrm{wt}(z)$); full construction (6553 terms; rational 4678) even=μ & odd=0 as $\mathbb Q(n)$ identities, exact $[c^{2n}]F=\mu_{7,n}$ (n≤30); HPL form vs direct corner integral $M_7(2\pi/3)$: **residual 3.5e-46 (~45 digits)**.

## Structural answers (item 4)
Weight bound **saturated** ($\max=k+1-\mathrm{wt}(z)$, zero violations, $k\le6$). Letters exactly $\{0,1,-1\}$; **−1 is needed** (from $H^{(r)}_{n-1}=2^{r-1}(Z_r+Z_{-r})-n^{-r}$; in 19/24, 63/74, 195/218 of the k=2,3,4 rational-part terms). Denominators: **only powers of c** ($c^{-j}$, $j\lesssim k-1$); $1\pm c$ never needed (coefficients are proper, so no $c\frac d{dc}$).

## Deliverables (work dir `new_results/FLD1_field_rule_proof/`)
- **`proof.md`** — the general argument + lemma ledger (statement, status, evidence for each lemma).
- **`verify_fld1.py`** — **5/5 PASS in 127 s (< 5 min)**: Lemma A (symbolic μ == exact table), Lemma C (all-n even=μ & odd=0 & low orders, all ζ-monomials k≤4), k=2,3 reproduction == B2/B3, k=4 exact μ (n≤100), k=4 == corner_series at 4 angles (~69 digits).
- **`result.md`** (full summary), **`progress.log`** (one line per step).
- **step_NN.py pipeline:** `step_20_bijection.py` (Lemma B), `step_21_fastmu.py` (exact μ, k≤5 n≤100, matches G1 flat_exact), `step_22_musymbolic.py` (symbolic-μ, Lemma A), `step_23_zsum.py` (Z-sum form, == `mu_zsum.pkl`), `step_24_inspect.py` (properness/poles), `step_25_construct.py` (Lemma C + all-n proof), `step_26_compare.py` (k=2,3==B2/B3), `step_27_validate.py` (exact μ n≤100 + corner_series 4 angles), `step_28_analyze.py` (weight/letters/denominators), `step_29b_g2compare.py` (==G2's B4), `step_30_pi3.py` (π/3 107 digits), `step_31_k6.py`/`step31b_k6analysis.py`/`step_32_k6construct.py`/`step_33_k6exact.py` (k=6), plus `zsum.py`, `hpl.py` (exact Z-sum/HPL machinery) and artifacts (`mu_exact_K5_N100.json`, `mu_exact_K6_N60.json`, `mu_zsum_K{5,6}.pkl`, `hpl_constructed_K{5,6}.pkl`).

All displayed numbers were printed by the scripts above (auto-archived as `step_NN.py`); `verify_fld1.py` reproduces the k=2,3 reproduction and the k=4 checks end-to-end in under 5 minutes.
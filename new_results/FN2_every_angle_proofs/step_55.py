result = r"""# FN2: Ahmed's Problem 3 at every angle — $\mathcal J_2,\mathcal J_3$ PROVED for general $f$; contact law PROVED

**Bottom line.** Both "derived and not proved" statements of the paper's Section sec:fn are
now **proved**, by transporting the $\theta=\pi/2$ proof of
`research/07_push/P3_fn_graph/claude/result.md` ([P3]) to every angle. The window
$\{y'>s\}$ does not see $\theta$; the only $\theta$-dependent inputs of that proof are
(i) the layer-cake density, (ii) the moments/cumulants of the log-slope, (iii) the tail
series of the density at $\infty$ — replaced here by Lemmas A, B, C. Full write-up in
**`proof.md`** (with a lemma ledger and two "$\theta$-audit" tables listing every replaced
step). Paper-ready statements in **`paper_proposal.tex`**. One-command re-verification
(<2 min) in **`verify_fn2.py`**.

## What is proved

- **Theorem 1 (every-angle expansion, general $f\in C^{K+2}$ near 1).**
  $2-L_n=\sum_{k\le K}\varepsilon^k\mathcal J_k(\Xi)+O(\varepsilon^{K+1}\log^{K+1}n)$,
  $\mathcal J_k=\sum_{i=1}^k\theta_{i,k-i}M_i(\Xi)+T_k$, with $M_i=2z\sum_r\binom ir
  m_r(\theta)\Xi^{i-r}$ and $T_k$ the slope-tail term with
  $d_l(c)=(-1)^l(1-c^2)P_{l-1}'(c)/l$. **Corollary:** $K=3$ gives exactly the paper's
  displayed $\mathcal J_1,\mathcal J_2,\mathcal J_3$ (sympy difference $0$ at general
  $(c,z)$; reduces to [P3]'s proved $Q_1,Q_2,Q_3$ at $\theta=\pi/2$; agrees with the exact
  series at $f=x$, also re-derived independently). So $\zeta(3)$ provably enters
  $\mathcal J_3$ at every angle, through $k_3(\theta)=2+2\cot^6\frac\theta2\log^3\cos^2
  \frac\theta2-12\cot^2\frac\theta2\,S_{1,2}(\sin^2\frac\theta2)$.
- **Theorem 2 (contact of order $m\ne1$ at every angle).** For
  $-\log f(1-\tau)=\beta\tau^m(1+O(\tau^\delta))$: (a) exact Lambert-$W$ reduction with
  relative error $O(n^{-\delta/m}\log^{1+\delta/m}n)$ (sharper $O(\cdot/\Lambda)$ version
  too); (b) $2-L_n=\frac{1-\cos\theta}{N}[F(\Lambda+k_1(\theta))+\frac{k_2(\theta)}2
  F''(\Lambda+k_1(\theta))+O(\Lambda^{1/m-3})](1+O(\delta_N))$;
  (c) the leading law and the **refined** law
  $$2-L_n=(1-\cos\theta)\Bigl(\frac{\log n}{m\beta n}\Bigr)^{1/m}
  \Bigl[1+\frac{(m-1)\log\log n+m\,k_1(\theta)+\log(m\beta)}{m\log n}
  +O\Bigl(\frac{\log\log n}{\log n}\Bigr)^2\Bigr],$$
  which explains the slow fall of the ratio to 1 ($\sim(\log\log n)/\log n$ from above for
  $m>1$).
- **Proposition ($u^\ast$ correction).** For the $m=2$ model integral,
  $(I-(1-c)u^\ast)u^\ast\to(1-c)k_1(\theta)/2$; the constant
  $\mathrm{Klayer}(c)=\frac12(\int_0^1 h_c/s-\int_1^\infty\rho_c/s)$ equals
  $(1-c)k_1(\theta)/2$ **exactly** (integration by parts). This corrects AH3's exploratory
  guess $\mathrm{Klayer}(c)-\frac{1-c^2}4$ (its step_29), which is off by $\frac{1-c^2}4$.

## Key new objects (proved)

- $\phi_\theta(s)=1-\frac{s+c}{\sqrt{1+2cs+s^2}}
  =\frac{1-c^2}{\sqrt{1+2cs+s^2}(\sqrt{1+2cs+s^2}+s+c)}$: decreasing, $0<\phi_\theta\le1-c$,
  mass $1-c$, $\phi_\theta=O(s^{-2})$ with explicit constant.
- Tail series $\phi_\theta(s)=\sum_{l\ge2}d_l(c)s^{-l}$ ($s>1$),
  $d_l(c)=(-1)^l(1-c^2)P_{l-1}'(c)/l$; at $c=0$ only even $l$, $d_{2m}(0)=-\binom{-1/2}m$
  (the $\theta=\pi/2$ coefficients).
- Moments: $\mathbb P(S_\theta\le s)\le s$, $\mathbb P(S_\theta>s)=O(s^{-1})$;
  $\mathbb E S_\theta^w<\infty$ iff $-1<w<1$; cumulants $k_j(\theta)=(-1)^jj!a_j$ from the
  paper's proved Legendre identity.

## Numerical evidence (all numbers printed by code)

| check | result | file |
|---|---|---|
| residuals after $\mathcal J_3$, scaled by $\varepsilon^4\log^4n$; 2 f × 2 θ (+π/2 control), $n=10^4,10^5,10^6$ | bounded: $|{\cdot}|\le0.022$ | `residuals_main.json`, `step_39` |
| resid after $\mathcal J_3$ $=\varepsilon^4\mathcal J_4(1+o(1))$ | ratios $0.99982\to0.999998$ | `residuals_K4.json`, `step_42` |
| contact ratios $D/L_0$ at $n=10^5$ | $1.1435,1.1598,1.1812$ (= AH3 exactly) | `contact_results.json`, `step_43` |
| contact ratios $D/L_1$ at $n=10^5$ | $0.9974,1.0116,1.0303$ (= AH3 exactly) | `step_43` |
| $D/$full-model at $n=10^4..10^8$ | $0.99991\to0.99999998$, err $\sim\log n/n$ | `step_43` |
| refined law (c): $(D/L_0-\mathrm{pred})(\Lambda/\log\Lambda)^2$ | in $[0.70,0.85]$ for $n=10^4..10^8$ | `step_44` |
| $u^\ast$ limit vs $(1-c)k_1/2$ vs AH3's guess | $\to(1-c)k_1/2$ (other off by $(1-c^2)/4$) | `step_45` |
| cumulants $k_{1,2,3}(\theta)$ vs closed forms | worst $6.1\cdot10^{-60}$, 7 angles | `step_29` |
| Mellin vs Legendre, $w=\pm0.95,\pm0.5,0.3$ | worst $7.5\cdot10^{-60}$, 4 angles | `step_29` |
| engine vs AH3 printed values | $\le2\cdot10^{-24}$ (their 23-digit rounding) | `step_36` |

## Labels

- **PROVED:** Lemma A, B, C, Theorem 1 + Corollary, Theorem 2, Proposition ($u^\ast$/Klayer).
  Algebra/identity steps are sympy-exact (difference $0$). The proof is by hand, written in
  full in `proof.md`; the only external inputs are [P3]'s proved propositions and the
  paper's proved Thm. P2/P2asym.
- **NUMERICAL:** the residual/ratio tables above (they confirm, but are not part of, the proofs).
- **OPEN:** nothing. (The uniformity of constants as $\theta\to0,\pi$ is not claimed; the
  theorems are for fixed $\theta\in(0,\pi)$.)

## Files

- `proof.md` — the deliverable: statements, proofs, $\theta$-audit tables, lemma ledger,
  numerical tables, proposed paper statements.
- `paper_proposal.tex` — the theorem + proposition for the paper (replaces both "derived"
  sentences and upgrades Table II).
- `verify_fn2.py` — consolidated verification, runs in ~100 s.
- `step_25`–`step_45` — the computational steps (see `README_STEPS.md`).
- `residuals_main.json`, `residuals_K4.json`, `contact_results.json` — data.
"""
with open("result.md","w") as f:
    f.write(result)

steps = r"""# Step files (auto-saved code blocks; run with the venv python, cwd=this folder)

Exploration/reading: step_01..step_23 (paper sec:fn extraction, P3 + AH3 + README reading,
setup of progress.log/result.md).

Analysis steps cited in proof.md:
- step_25: Lemma B symbolic — phi_theta alt form, monotonicity, tail d_l(c) = (-1)^l(1-c^2)P'_{l-1}(c)/l,
  c=0 -> -C(-1/2,m)   [sympy, exact]
- step_26: Lemma B numeric — mass 1-c, explicit s^-2 bound, tail remainder constants
- step_29: Lemma C stable — cumulants k1,k2,k3(theta) vs closed forms (worst 6.1e-60, 7 angles);
  Mellin vs Legendre (worst 7.5e-60); k4(theta) table
- step_30: Theta(W,mu) series inversion, theta_ij (i+j<=3), theta_{1j} formula
- step_31: MASTER FORMULA -> paper's displayed J1,J2,J3 (sympy difference 0 at general (c,z))
- step_32: reductions: theta=pi/2 -> P3's Q1,Q2,Q3; f=x -> exact_J's (difference 0)
- step_33: independent f=x every-angle series from layer cake + CGF + tail (difference 0)
- step_36: D_n engine validation vs AH3 printed values (<=2e-24)
- step_39: residuals after J3 scaled by eps^4 log^4 n, 2 f x 3 th x n=1e4..1e6 (bounded) [data: residuals_main.json]
- step_40: theta_{4j} symbolic, matches P3
- step_41: T4 general (c,z); at c=0,z=1/2 matches P3's T4
- step_42: K=4 numerics: resid3/(eps^4 J4) -> 1; resid4/(eps^5 log^5 n) bounded [data: residuals_K4.json]
- step_43: contact numerics: sin(pi x/2), 3 angles, n=1e4..1e8; D vs L0/L1/L2/L3/full model [data: contact_results.json]
- step_44: refined leading law quantitative; K(c) numeric (first pass)
- step_45: u* limit: (I-(1-c)u*)u* -> (1-c) k1/2 (not Klayer-(1-c^2)/4)

Consolidated re-run of everything: verify_fn2.py (~100 s).
"""
with open("README_STEPS.md","w") as f:
    f.write(steps)

import time
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} proof.md complete (Theorems 1-2, Corollary, Proposition, ledger, paper proposal); result.md final; README_STEPS.md; verify_fn2.py OK in 102s  [DONE]\n")
print("result.md, README_STEPS.md written; progress.log updated")
import os
print(sorted(os.listdir(".")))

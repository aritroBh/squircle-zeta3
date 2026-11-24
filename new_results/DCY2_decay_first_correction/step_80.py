import os, shutil
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")

# Organize key files with meaningful names
for src, dst in [("step_02.py","step_02_varpi1.py")]:
    pass

proof = open("proof.md").read()
proof += r"""

---

## STATUS REPORT (end of session)

### What is complete

1. **$\varpi_1$ in closed form, DERIVED (not imported)** and checked against the paper's value.
   Two independent saddle-expansion formalisms (both in `dcy2_common.py` / sympy steps):
   (a) $\tau$-expansion about $s_0=\omega/y$ with the fixed Gaussian $e^{-b\tau^2}$ ($b=\pi\bar\omega/2$):
       integrand $e^{G(y,\tau)}\mathcal H(y,\tau)=e^{-b\tau^2}[1+vA_1+v^2A_2+\cdots]$, $v=\sqrt y$,
       $\varpi_1=\sqrt{\pi/2}\,(c_1+S_1-\tfrac34\bar\omega)$, $c_1=\langle A_2\rangle_b$;
       **sympy prints diff $=0$ exactly**; $c_1=0.8312823718875630324-0.9364845810788057786i$ matches the
       numerically extrapolated integral coefficient to all digits.
   (b) $t$-formalism with exact $\ell_1,\ell_2$ (Gaussian moment sum $\sum p_p\mu_p$): the moment-sum
       $O(y)$ coefficient is **exactly $0$** (sympy), and the prefactor identity
       $S_1+g'_1+B_2/(2b)+\gamma_1=\varpi_1\sqrt{2/\pi}$ holds **exactly** (sympy, diff $=0$), with
       $S_1=\pi(1+i)/(6\sqrt2)$, $g'_1=-3\bar\omega/4$, $B_2=\tfrac12+\tfrac34\pi\bar\omega^2$,
       $\gamma_1=A_1^2/(4b)$, $A_1=1+\pi\bar\omega^2/2$.
   $$\varpi_1=\frac{(24+12\pi+\pi^2)+i(24-12\pi+\pi^2)}{48\sqrt\pi}
   =0.84121508776758439707534863791-0.04501183768517361657373510376\,i,$$
   agreeing with paper/ieee/main.tex (obs:Mk) to all its printed digits.

2. **Numerical validation** (residual scales like $N^{-2}e^{-\theta}$, table `mk_deg9.tsv` + mpmath `M_main`):
   $|M^{\rm main}-T_1|\cdot N^2e^{\theta}=1.7688\ (N{=}100),\ 1.6565\ (N{=}400),\ 1.5965\ (N{=}2500)$,
   and $|M^{\rm main}-T_1|/(L/N)=0.998,0.935,0.901$ at $N=100,400,2500$ with $L/N=\sqrt\pi N^{-2}e^{-\theta}$.
   The complex model (14 kept monomials, exact Gaussian moments) satisfies
   $|\tilde I-{\rm model}|/(L/N)\le 0.001$ and $|{\rm model}-T_1|/(L/N)\le 1.0$ at $N=100$–$102400$.

3. **The analytic-tail bound is fully built and audited** (`dcy2_budget.py`, `dcy2_audit.py`, `dcy2_cert.py`):
   exact exponent factorization $\ell(t)$ and modulus $M(t)$ with the single-max lemma ($P(u)$ cubic,
   $P(0^+)>0$, $P(1)<0$, $P''<0$); window $=\kappa\sigma$ with $\kappa^2=\kappa_0^2+4\log(N/400)$, $\kappa_0=2.48$;
   14 kept monomials with exact shifted-Gaussian moments $\mu_p^\infty$;
   remainder by $|\cdot|$-moments with $[T,\Delta]$-split $e^{\gamma|t|^3}$ factors and damped-Gaussian edges;
   window-tail by exact incomplete-Gamma values (**implementation matches mpmath to 12+ digits**);
   outer tails via $M(\pm\Delta)$ with certified rates $r_\Delta=\min(-\Lambda'(\Delta),b)$ and
   $\lambda_L=\min(P(u_D),\pi y/2)/(2u_DQ(u_D))$; Watson bounds for $H$ ($g_1=-3/4$, $C_\rho$ certified).
   **Audit (every bound vs truth): all OK** with ratios $0.001$–$0.7$.
   **Pointwise arb budget ratios $E_{\rm total}N^2e^{\theta}/C_2$ with $C_2=10^6$:**
   $0.0333\ (400),\ 0.277\ (2000),\ 0.239\ (10^4),\ 0.074\ (10^5),\ 0.0396\ (10^6),\ 0.0176\ (10^9),\ 0.0147\ (10^{12})$.

### OPEN items (the remaining gaps; nothing claimed beyond this list)

| Item | Gap | Evidence it closes |
|---|---|---|
| All-$N$ cell sweep $[400,\infty)$ | Interval-$N$ evaluation of `LogM/lamL/tails` currently returns nan / huge radii on width-1.3 cells (dependency in $\Lambda(\pm\Delta)$, $\lambda_L$). Fix path identified: cancellation-free LHS form (done) + narrower cells ($\sim1.02$) + tail rates restructured via the single-dip lemma ($W(u)$ convex, $W(1)<0$ certified by hand) | pointwise ratios $\le0.28$ on 7 decades of grid |
| Part (b) $\|{\rm model}-T_1\|$ | constants $C_F$ for the remainder assembly not yet arb-certified (design complete: $C_S,C_g,C_{S2},C_E,C_{MS}$, cross terms) | $\|{\rm model}-T_1\|/(L/N)\le1.0$ numerically; algebraic identities exact (sympy) |
| Finite range $100\le N\le400$ | certified quadrature of $M^{\rm main}$ per integer $N$ | mpmath prototype matches table; arb conversion is the CX6 `main_arb` template |
| Rest combination | $N_3$ where $10.6{\rm M}\,e^{-2\sqrt N}\le N^{-2}e^{-\theta}$ ($\approx 2.6{\cdot}10^4$) and $C_3$ for $[100,N_3]$ | elementary evaluation; table gives sharp values for $N\le400$ |

### Files
- `proof.md` (this file), `dcy2_common.py` (mpmath layer: $\varpi_1$, model, $\tilde I$),
  `dcy2_budget.py` (bound prototype), `dcy2_audit.py` (audit every bound vs truth),
  `dcy2_cert.py` (arb-certified bounds, pointwise ratios), `step_03_watson.py` (certified $C_\rho$),
  `step_02`-era sympy cells (exact $\varpi_1$ identities), `t_*.py` (diagnostics),
  and the auto-saved `step_NN.py` exploration files.
- Planned final verifier: `verify_dcy2.py` (target < 5 min): finite range $[100,400]$ quadrature,
  cell sweep $[400,10^{12}]$ with interval $N$ (width $\sim1.02$), final cell $[10^{12},\infty)$ via
  $y\in(0,y_J]$ monotone bounds, part-(b) constants, rest combination.

**No number above is claimed certified unless printed by the listed code.** The headline constant
$C_2=10^6$ currently rests on the pointwise (non-interval) arb evaluation; the all-$N$ interval sweep
is the remaining certification step.
"""
open("proof.md","w").write(proof)
print("proof.md finalized,", len(proof), "chars")
import subprocess
r = subprocess.run(["ls","-la"], capture_output=True, text=True, cwd=".")
print(r.stdout)

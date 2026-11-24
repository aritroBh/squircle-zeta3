WD='<project>/research/08_theory/FLD1_field_rule_proof'
s=open(WD+'/proof.md').read()
i=s.find('## Lemma D'); j=s.find('## Lemma E')
old = s[i:j]
new = """## Lemma D (weight bound)

$\\operatorname{weight}(z_j)+|w_j|\\le k+1$ for every term of the constructed form, for **all $k$**.

**Proof (repaired pole-budget, after the referee).** For a nested-sum term $\\rho(n)Z_w(2n)$ and a point $p$, let
$\\operatorname{ord}_p(\\rho)$ be the pole order of $\\rho$ at $p$ (0 if regular), and let $\\nu(\\rho)=\\max_p\\operatorname{ord}_p(\\rho)$.
We show, at *every* location $p$,
$$|w|+\\operatorname{ord}_p(\\rho_w)\\ \\le\\ k+1-\\operatorname{weight}(z)\\qquad(\\ast)$$
for every term $\\rho_w(n)Z_w(2n)$ of $[\\mu_{k,n}]_z$; then $|w|+\\nu(\\rho_w)\\le k+1-\\operatorname{weight}(z)$.

The sources of poles in $\\rho_w$ (Lemma A) split into two disjoint structural classes and the polygamma shifts:
* **the prefactor** $-1/(n(2n-1))$: simple poles at $\\{0,\\frac12\\}$ only;
* **the Pochhammer ratio** $(n)_j/(2n)_j=\\prod_{t=0}^{j-1}\\frac{n+t}{2n+t}$: simple poles at $n=-t/2$, $t\\ge1$ only —
  crucially, *no pole at $n=0$*, because the $t=0$ factor is $\\frac{n}{2n}=\\frac12$ (removable). Its pole set
  $\\{-t/2:t\\ge1\\}$ is therefore **disjoint** from the prefactor's $\\{0,\\frac12\\}$;
* **the polygamma shifts** $(n+t)^{-r}$ and $(2n+t)^{-r}$ from Lemma A(iv): pole order $r$ at $n=-t$ (integer) or
  $n=-t/2$ (half-integer). A shift arises exactly when a polygamma $\\psi^{(r-1)}$ is *not* kept as a full harmonic
  $H^{(r)}\\mapsto Z_{\\pm r}$; kept, it contributes $r$ to $|w|$ and $0$ to the pole order; shifted, it contributes $0$
  to $|w|$ and $r$ to the pole order. Hence the total shift order, at any single location, is bounded by the part of the
  harmonic budget *not* used in $|w|$, namely $\\le (k-\\operatorname{weight}(z))-|w|$ (Lemma A1; the budget is additive
  under stuffle, which preserves $|w|$ and adds pole orders).

Now at each location:
* at $p\\in\\{0,\\frac12\\}$: $\\operatorname{ord}_p\\le 1+\\big[(k-\\mathrm{wt}(z))-|w|\\big]$ (prefactor's single pole + shift budget);
* at $p=-t/2$ ($t\\ge1$): $\\operatorname{ord}_p\\le 1+\\big[(k-\\mathrm{wt}(z))-|w|\\big]$ (Pochhammer's single pole + shift budget);
* at $p=-t$ ($t\\ge1$ integer): $\\operatorname{ord}_p\\le (k-\\mathrm{wt}(z))-|w|$ (shift budget only; no structural pole).
In all three cases $|w|+\\operatorname{ord}_p(\\rho_w)\\le |w|+1+(k-\\mathrm{wt}(z)-|w|)=k+1-\\mathrm{wt}(z)$, proving $(\\ast)$.
(The "$+1$" is incurred at most once per location, since prefactor and Pochhammer poles are disjoint.)

Finally, in Lemma C the HPL word for a partial-fraction term of pole order $\\kappa$ is $w_0=\\mathrm{bji}(1,\\kappa,w)$ of
weight $|w_0|=\\kappa+|w|\\le\\nu(\\rho_w)+|w|$ (Lemma B2), so
$$\\operatorname{weight}(z)+|w_0|\\ \\le\\ \\operatorname{weight}(z)+\\big(k+1-\\operatorname{weight}(z)\\big)=k+1.$$
And the $c\\frac d{dc}$ branch never raises the word length: $c\\frac d{dc}H_w=c f_{w_1}H_{w'}$ equals $H_{w'}$
(weight $|w|-1$) if $w_1=0$, and $(\\frac1{1\\mp c}-1)H_{w'}$ — a prefix-sum HPL of weight $|w'|+1=|w|$ — if $w_1=\\pm1$.
So the bound is unaffected by the (vacuous for $k\\le7$) polynomial-part clause. $\\blacksquare$

**Status.** PROVED for all $k$ (the repair accounts for the Pochhammer poles, which FLD1's budget had forgotten).
Verified per word *and per pole location* at $k=7$: 0 violations, saturated (`FLD1R step_34`; sharpening our
`step_28`/budget check for $k\\le6$).

---

"""
s = s[:i] + new + s[j:]
open(WD+'/proof.md','w').write(s)
print("Lemma D rewritten for all k")
WD='<project>/research/08_theory/FLD1_field_rule_proof'
s=open(WD+'/proof.md').read()
i=s.find('## Lemma ledger')
newtail = """## Lemma ledger

| Lemma | Statement | Status | Evidence |
|---|---|---|---|
| thm:flat | $M_k=\\sum\\mu_{k,n}c^{2n}$, $\\mu_{k,n}\\in\\mathbb Q[\\zeta(2..k)]$ | PROVED (paper) | cited |
| A (formula, A1) | $\\mu_{k,n}$ = nested sums; $|w|\\le k-\\mathrm{wt}(z)$; $i+\\ell=k$ exactly | PROVED (all $k$; verified $k\\le7$) | `step_21/22/23`; referee `step_26/28/31/42` |
| A2 (poles) | poles in $\\{\\frac12,0\\}\\cup\\{-\\frac t2:t\\le2k-2\\}\\subset\\frac12\\mathbb Z$ | PROVED (coarser form, all $k$); sharp set $\\{t\\le k-1\\}$ PROVED FOR $k\\le7$ | `step_24/31`; referee `step_32/34` |
| A2 (properness) | $\\deg\\rho_w=-2$ (denominators powers of $c$) | PROVED FOR $k\\le7$; OPEN for $k\\ge8$ (needs $Z$-sum independence lemma) | `step_24/31`; referee `step_32/34` (553 word-terms, all deg $-2$) |
| B | bijection $w\\leftrightarrow(\\varepsilon,k,v)$; $|w|=k+|v|$ | PROVED (all $k$) | `step_20` (242 words); referee re-verified independently (728 words len≤6, `step_24`) |
| C | even nested sums $\\to$ HPL (algorithm terminates, is correct; regularity at $n=1$) | PROVED (all $k$); verified by construction FOR $k\\le7$ | `step_25/32`; referee `step_35/36` (k=7) |
| D | $\\mathrm{wt}(z)+|w|\\le k+1$ | **PROVED for all $k$** (repaired pole-budget) | this proof; referee `step_34` (per word & pole, 0 violations, saturated) |
| E | letters $\\{0,\\pm1\\}$; $-1$ needed | PROVED (all $k$) | `step_28`; at $k=7$, $-1$ in $4513/4678=96.5\\%$ of rational terms (referee) |
| F | $S$-units $\\langle2,c,1-c,1+c\\rangle$ | PROVED | symbol letters; paper's $\\mathrm{Li}_3$ arguments; referee (MPL args $\\pm c,\\pm1$) |

**The theorem (existence of the HPL form, letters, $S$-units, weight bound) is PROVED for all $k$.** Lemmas A, B, C, E, F
are written for general $k$; Lemma C's polynomial-part clause makes existence independent of the (unproven for $k\\ge8$)
properness, costing at most $(1\\pm c)$ denominators; and after the referee's repair, Lemma D's weight bound is itself
proved for all $k$. The only clauses not proved for all $k$ are the **power-of-$c$ denominators** and the **sharp pole
set**, both PROVED FOR $k\\le7$ and open for $k\\ge8$ (they would follow from properness via a $Z$-sum independence lemma).
The whole construction (symbolic $\\mu$, $Z$-sum form, HPL form, all-$n$ identity) is verified exactly for $k\\le7$; an
independent referee (FLD1R) re-verified Lemmas A, B, C, E, F for all $k$ and supplied the repairs to A2, C and D
incorporated here.

---

## Verification (see verify_fld1.py and result.md)

* **$k=2,3$ reproduction.** The construction reproduces Theorem thm:M2alpha (B2, 24 terms) and B3 (74/15/5 terms)
  *exactly* (Taylor agreement to $c^{40}$; `step_26`), and the all-$n$ Z-sum identities (even part $=\\mu$, odd part $=0$)
  hold (`step_25`).
* **$k=4,5$ construction + validation.** HPL forms built by the algorithm: exact $[c^{2n}]=\\mu_{k,n}$ for $n\\le100$ (all
  $\\zeta$-monomials, odd and non-positive parts $0$; `step_27`), and agreement with the corner integral
  (`corner_series`) at 4 angles to **$\\ge40$ digits** ($\\ge76$ digits at $2\\pi/3,\\pi/2,3\\pi/4$ at N=600/dps=80; and
  **107 digits** at $\\pi/3$ with N=2500/dps=110, confirming the form is exact there; `step_27`, `step_30`).
* **Consistency with G2_M4alpha.** Our $k=4$ equals G2's independently-fitted-and-proved B4 (R4, $\\Psi_2$,
  $-\\zeta(3)\\Psi$, $\\zeta(4)F_1$, $\\zeta(2)^2F_2$) as functions (`step_29`, `step_29b`).
* **$k=6$ (beyond the required range).** Symbolic $\\mu_{6,n}$ == exact table ($n\\le60$); Z-sum form proper (deg $-2$),
  weight budget $|w|+\\nu\\le 7-\\mathrm{wt}(z)$ (saturated), poles $\\{-\\frac52,\\dots,0,\\frac12\\}$; and the full
  construction + all-$n$ identity (even $=\\mu$, odd $=0$) for all 11 $\\zeta$-monomials (`step_31`, `step_32`).
* **$k=7$ (referee FLD1R).** Symbolic $\\mu_{7,n}$ (429 monomials) == exact table ($n\\le30$, two independent
  implementations identical); the Z-sum canonical form (553 word-terms) is proper (all $\\deg=-2$) with poles
  $\\{-3,-\\frac52,-2,-\\frac32,-1,-\\frac12,0,\\frac12\\}$; weight budget per word *and per pole location* has 0 violations
  (saturated at $8-\\mathrm{wt}(z)$); the full construction (6553 HPL terms; rational part 4678 terms) satisfies even $=\\mu$
  and odd $=0$ as $\\mathbb Q(n)$ identities and exact $[c^{2n}]F=\\mu_{7,n}$ for $n\\le30$; and the HPL form agrees with the
  direct corner integral $M_7(2\\pi/3)$ to residual $\\mathbf{3.5\\cdot10^{-46}}$ ($\\sim45$ digits, dps=50, N=500)
  (referee `step_26`–`step_39`). This is the computational content of "power-of-$c$ denominators PROVED FOR $k\\le7$".
"""
s = s[:i] + newtail
open(WD+'/proof.md','w').write(s)
print("ledger + verification updated")
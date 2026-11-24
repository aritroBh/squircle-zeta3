import os
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")

md = r"""# OCT5D2 — the depth test for the regular octagon's P⁽⁵⁾(1)

**Headline (symbol level, PROVED; number level modulo the period conjecture):**
the weight-5 part of the regular octagon's \(P^{(5)}(1)/u\) has a **nonzero** antisymmetrised (3,2)
cobracket component. The regular octagon **does need a depth-two polylogarithm at order five** —
the same phenomenon as the squircle's \(M_4\). Contrary to the OCT5D planner's heuristic, the
obstruction does **not** come from the mixed two-ratio block alone: the single-ratio block also
carries depth two. The obstruction pair is **not proportional** to \(M_4\)'s \((84,-38\sqrt2)\), so the
octagon's depth-two correction is not a scalar multiple of \(M_4\)'s \(N_5\).

## Setup (reused, validated)

- Exact edge integrand: 461 monomials in \((\log t,\log(1-t),\log(1+\kappa_0 t),\log(1+\kappa_1 t),\log2,\log(1+\sqrt2))\),
  coefficients \(u\cdot\mathbb Q(\sqrt2)\), from `OCT5D_octagon_order5_depth/parsed_mons.pkl`;
  \(P^{(5)}(1)/u = 960\sum \text{coeff}\cdot(\log2)^e(\log(1+\sqrt2))^f\int_0^1 t^n(\log t)^a(\log(1-t))^b(\log(1+\kappa_0t))^c(\log(1+\kappa_1t))^d dt\).
- Weight-5 part = monomials with \(a+b+c+d+e+f=5\); classes: local (c=d=0, 56 monomials),
  single-ratio (exactly one of c,d nonzero, 140), mixed (c,d≥1, 56).
- Letters for the \(G\)-rewrite: \(0,1,p=2+\sqrt2=-1/\kappa_0,\ q=-(1+\sqrt2)=-1/\kappa_1\).
  All letter differences are S-units in \(\{2,1+\sqrt2\}\), so the full Goncharov symbol lives in
  \(V^{\otimes 5}\), \(V=\langle X=\log2,\ Y=\log(1+\sqrt2)\rangle\) (2-dimensional, as for \(M_4\)).

## Task 1 — Rewrite (PROVED algebra / NUMERICAL check)

Each integral is rewritten **exactly** over \(\mathbb Q(\sqrt2)\) as
\(\int_0^1 t^n(\log t)^a(\log(1-t))^b(\log(1+\kappa_0t))^c(\log(1+\kappa_1t))^d dt
= a!b!c!d!\sum_{w\in 0^aш1^bш p^cш q^d} I_G(n,w)\),
\(I_G(n,w) = \frac{1-w_1^{n+1}}{n+1}G(w;1) - \frac1{n+1}\sum_{j=0}^n w_1^{\,n-j}I_G(j,w[1:])\), \(I_G(n,())=1/(n+1)\),
with shuffle-regularised \(G_*\) (trailing-zero words; words starting with 1 provably never get a
nonzero coefficient). 756 integrals, 1024 distinct words, weights ≤5 (`rewrite_exact.pkl`).

Numerical validation of the rewrite (G-table: two-centre power series at \(t=1/2\) + Chen
concatenation, dps 170, truncation N=480; independently validated against `mpmath.polylog`,
direct quadrature, nested quadrature, shuffle identities — residuals \(10^{-147}\)…\(10^{-171}\)):

| block | rewrite value minus itab110-assembled target | digits |
|---|---:|---:|
| local | \(2.8\times10^{-102}\) | ≥102 |
| single-ratio | \(3.3\times10^{-106}\) | ≥106 |
| mixed | \(2.6\times10^{-107}\) | ≥107 |
| mixed vs task-quoted 4.2299562389856442293498943256846704703130475386317 | \(-2.6\times10^{-50}\) | ≥49 (quoted value truncated there) |

A second, independent **dual rewrite** (substitution \(u=1-t\), mapping
\(\log(1+\kappa_0t)\mapsto -\tfrac12\log2+G(q;u)\), \(\log(1+\kappa_1t)\mapsto\tfrac12\log2+G(p;u)\)) reproduces
the same integrals numerically (worst \(1.2\times10^{-80}\) on 40 random integrals at dps 80) and the
**same top-weight symbols exactly for all 756 integrals**.

## Task 2 — Symbol and controls (PROVED)

Full number-level Goncharov symbol with every constant entry \((a_i-a_j)\), exact S-unit
factorisation over \(\mathbb Q(\sqrt2)\) modulo torsion (the validated `g_symbol` of
`A3_W4_proof/w4_symbol.py`; same machinery and conventions as M4c). Controls:

- **Li₅ at 7 S-units** (\(b,-b,-b^2,\sqrt2/2,1/2,1-\sqrt2/2,2-\sqrt2\)) pairs to \((0,0)\): PASS.
- **300 random shuffle products** of weights (1,4),(2,3) pair to \((0,0)\): PASS.
- **M4c's \(W_5\)** (from `symbol.json`) pairs to **\((84,-38\sqrt2)\)**: PASS (both with
  `w5_symbol.obstruction` and an independent planner-style re-pairing).
- **G-symbol identities**: \(S\,G(0^{k-1},p;1) = -S\,\mathrm{Li}_k(1/p)\) exactly for \(k=2..5\);
  shuffle identity \(S\,G(p;1)ш S\,G(q;1) = S\,G(p,q;1)+S\,G(q,p;1)\): PASS (the underlying
  numerical values also checked against quadrature/polylog at \(10^{-147}\)…\(10^{-171}\)).
- **Scaling invariance** \(S\,G_0(w;1) = S\,G_0(\lambda w;\lambda)\) for 780 convergent words
  (\(\lambda = 1/p\) and \(\lambda=-1/q\), all letters \(\{0,1,p,q\}\)): PASS, 0 mismatches.
  (Regularised trailing-zero words are convention-fixed, same convention as M4c; the rewrite's
  numerical check validates the convention end-to-end.)
- **N₅ reproduced**: the four \(M_4\) atoms \(\mathrm{Li}_{3,2}(b,1),\mathrm{Li}_{3,2}(1/\sqrt2,1),
  \mathrm{Li}_{3,2}(-b,1),\mathrm{Li}_{3,2}(b,-1)\) with cobracket coefficients pair to
  \((84,-38\sqrt2)\): PASS — independent confirmation that the pairing functional is calibrated.

## Task 3 — The depth test (PROVED, symbol level)

Pairings against \([X,[X,Y]]\otimes[X,Y]\) (first) and \([Y,[X,Y]]\otimes[X,Y]\) (second),
exactly as M4c. Results (elements of \(\mathbb Q(\sqrt2)\)):

| block | pair 1 | pair 2 |
|---|---|---|
| local | \(0\) | \(0\) |
| single-p (c≥1,d=0) | \(\tfrac{15}{8}(-3597+2544\sqrt2)\) | \(\tfrac{45}{4}(-454+321\sqrt2)\) |
| single-q (c=0,d≥1) | \(\tfrac{15}{8}(-3597+2544\sqrt2)\) | \(\tfrac{45}{4}(-454+321\sqrt2)\) |
| single (total) | \(\tfrac{15}{4}(-3597+2544\sqrt2)\) | \(\tfrac{45}{2}(-454+321\sqrt2)\) |
| mixed | \(\tfrac{45}{4}(1067-189\sqrt2)\) | \(\tfrac{45}{2}(494-345\sqrt2)\) |
| **total W₅/u** | \(\mathbf{-1485+1035\sqrt2}\) | \(\mathbf{900-540\sqrt2}\) |

(single-p and single-q are exactly equal — the \(\kappa_0\leftrightarrow\kappa_1\), \(u=1-t\) symmetry.)

**The total obstruction \((-1485+1035\sqrt2,\ 900-540\sqrt2)\) is nonzero** in both components:
the octagon's \(P^{(5)}(1)\) needs a depth-two polylogarithm at order five (symbol level; number
level modulo the period conjecture). Two independent assemblies (direct and \(u=1-t\)-dual) give
identical symbols and obstructions.

**Ratio to \(M_4\)'s \((84,-38\sqrt2)\):** the pair is **not proportional**. Component ratios in
\(\mathbb Q(\sqrt2)\):
\[
\frac{-1485+1035\sqrt2}{84}=\frac{15(-33+23\sqrt2)}{28},\qquad
\frac{900-540\sqrt2}{-38\sqrt2}=\frac{45(6-5\sqrt2)}{19},
\]
which are unequal (numerically \(-0.248\ldots\) vs \(-2.537\ldots\)). So the octagon's depth-two
correction is a **different combination** of depth-two atoms, not \(\text{ratio}\cdot N_5\).

## Task 4 — Depth-two correction and numerical form (in progress)

The \(M_4\) atom set \(\mathrm{Li}_{3,2}(x,y)\), \(x\in\{b,1/\sqrt2,-b\}\), \(y=\pm1\), has
(3,2) pairings spanning the 2-dimensional obstruction space, so a pair-matching correction exists
within the set (e.g. \(8(1485-1035\sqrt2)\,\mathrm{Li}_{3,2}(1/\sqrt2,1) +
(-3600+2160\sqrt2)\,\mathrm{Li}_{3,2}(b,1)\) matches the pair). The pair alone does not fix the
correction uniquely (pair-kernel); the numerical classical-form test over a combined
\{classical\} \(\cup\) \{\(\mathrm{Li}_{3,2}\) atoms\} basis at high precision is the decisive check.
Status: running.

## Labels

- PROVED (exact symbol algebra): the rewrite algebra, all controls, the obstruction pairs,
  non-proportionality.
- NUMERICAL (digits): the rewrite-vs-value checks (\(\ge 49\)…\(\ge 107\) digits), the \(G\)-table
  validations, any closed-form statement in Task 4.
"""
open("result.md","w").write(md)
print("result.md written, %d chars" % len(md))

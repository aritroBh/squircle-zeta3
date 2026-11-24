# OCT5D2 — Depth test for the regular octagon's P⁽⁵⁾(1): RESULT

## Bottom line

**The regular octagon needs a depth-two polylogarithm at order five** (PROVED at symbol level; number level modulo the period conjecture). The antisymmetrised (3,2) pairings of the weight-5 part of \(P^{(5)}(1)/u\), \(u=\sqrt{2+\sqrt2}\), against \([X,[X,Y]]\otimes[X,Y]\) and \([Y,[X,Y]]\otimes[X,Y]\) (X = log 2, Y = log(1+√2)) are

\[
\boxed{\;\bigl(-1485+1035\sqrt2\bigr)\;=\;45(-33+23\sqrt2),\qquad \bigl(900-540\sqrt2\bigr)\;=\;180(5-3\sqrt2)\;}
\]

both **nonzero**. The pair is **not proportional** to M₄'s \((84,-38\sqrt2)\): component ratios \(15(-33+23\sqrt2)/28\) and \(45(6-5\sqrt2)/19\) in \(\mathbb Q(\sqrt2)\) are unequal (numerically \(-0.248\) vs \(-2.537\)). Two independent exact assemblies (direct rewrite and \(u=1-t\) dual rewrite) give identical symbols and pairings. Contrary to the OCT5D planner's heuristic ("the mixed two-ratio block is the only block that can carry depth two"), **the single-ratio block also carries depth two**, and it does so in two exactly equal halves (\(\kappa_0\leftrightarrow\kappa_1\) symmetry):

| block | pairing 1 (vs \([X,[X,Y]]\otimes[X,Y]\)) | pairing 2 (vs \([Y,[X,Y]]\otimes[X,Y]\)) |
|---|---|---|
| local (c=d=0) | \(0\) | \(0\) |
| single-p (c≥1,d=0) | \(\tfrac{15}{8}(-3597+2544\sqrt2)\) | \(\tfrac{45}{4}(-454+321\sqrt2)\) |
| single-q (c=0,d≥1) | \(\tfrac{15}{8}(-3597+2544\sqrt2)\) | \(\tfrac{45}{4}(-454+321\sqrt2)\) |
| mixed (c,d≥1) | \(\tfrac{45}{4}(1067-189\sqrt2)\) | \(\tfrac{45}{2}(494-345\sqrt2)\) |
| **total W₅/u** | \(\mathbf{-1485+1035\sqrt2}\) | \(\mathbf{900-540\sqrt2}\) |

## Task 1 — Rewrite (PROVED algebra; NUMERICAL check ≥ 49–107 digits)

Each of the 756 weight-5 monomial integrals is rewritten **exactly** over \(\mathbb Q(\sqrt2)\):
\(\int_0^1 t^n(\log t)^a(\log(1-t))^b(\log(1+\kappa_0t))^c(\log(1+\kappa_1t))^d dt = a!\,b!\,c!\,d!\sum_{w\in 0^aш 1^bш p^cш q^d} I_G(n,w)\) with \(p=2+\sqrt2=-1/\kappa_0\), \(q=-(1+\sqrt2)=-1/\kappa_1\), \(I_G(n,w)=\frac{1-w_1^{n+1}}{n+1}G(w;1)-\frac1{n+1}\sum_{j=0}^n w_1^{n-j}I_G(j,w[1:])\), \(I_G(n,())=1/(n+1)\), shuffle-regularised \(G_*\) (words starting with 1 provably get zero coefficient). 1024 distinct words, all weights ≤ 5. Numerical G-table (two-centre power series at \(t=1/2\) + Chen concatenation, validated independently against `mpmath.polylog`, direct and nested quadrature, shuffle identities: residuals \(10^{-147}\)–\(10^{-432}\)):

| check | residual | digits |
|---|---:|---:|
| mixed block vs itab110-assembled 110-digit value 4.229956238985644229349894325684670470313047538631673924… | \(2.6\times10^{-107}\) | 107 |
| mixed block vs task-quoted 4.2299562389856442293498943256846704703130475386317 | \(-2.6\times10^{-50}\) | 49 (quoted truncation) |
| single-ratio block vs 110-digit target | \(3.3\times10^{-106}\) | 106 |
| local block vs 110-digit target | \(2.8\times10^{-102}\) | 102 |
| all 756 individual integrals vs itab110 | worst \(2.9\times10^{-86}\) (at dps 90), \(2.0\times10^{-102}\) (at dps 170) | ≥86 |
| dual rewrite (\(u=1-t\)) vs itab110, 40 random integrals | \(1.2\times10^{-80}\) | 80 |

## Task 2 — Symbol and controls (PROVED)

Full number-level Goncharov symbol with **every** constant entry \((a_i-a_j)\), exact S-unit factorisation over \(\mathbb Q(\sqrt2)\) modulo torsion (validated `g_symbol` from `A3_W4_proof/w4_symbol.py`; identical machinery and conventions as M4c). **All controls PASS:**
- Li₅ at 7 S-units (\(b,-b,-b^2,\sqrt2/2,1/2,1-\sqrt2/2,2-\sqrt2\)) pairs to \((0,0)\);
- 300 random shuffle products pair to \((0,0)\);
- M4c's \(W_5\) (from `symbol.json`) pairs to \((84,-38\sqrt2)\), with both `w5_symbol.obstruction` and an independent planner-style re-pairing;
- G-symbol test-word identities: \(S\,G(0^{k-1},p;1)=-S\,\mathrm{Li}_k(1/p)\) exactly for \(k=2..5\), \(S\,G(p;1)ш S\,G(q;1)=S\,G(p,q;1)+S\,G(q,p;1)\) (numerical values also quadrature-checked);
- scaling invariance \(S\,G_0(w;1)=S\,G_0(\lambda w;\lambda)\) for 780 convergent words (\(\lambda=1/p,-1/q\)): 0 mismatches (regularised trailing-zero words are convention-fixed, same convention as M4c, validated end-to-end by the rewrite's numeric check);
- **M4c's \(N_5\) reproduced**: the four atoms \(\mathrm{Li}_{3,2}(b,1),\mathrm{Li}_{3,2}(1/\sqrt2,1),\mathrm{Li}_{3,2}(-b,1),\mathrm{Li}_{3,2}(b,-1)\) with cobracket coefficients pair to \((84,-38\sqrt2)\);
- direct and \(u=1-t\)-dual rewrites give identical top-weight symbols for **all 756 integrals**.

## Task 3 — Depth test (PROVED symbol level)

As tabulated above: the total obstruction is nonzero in both components ⇒ **depth two is needed at order five**, as for the squircle's \(M_4\). The obstruction is not a scalar multiple of M₄'s \((84,-38\sqrt2)\).

## Task 4 — Correction and numerical status

- **PROVED (exact symbol algebra):** \(S(W_5^{oct})\) lies in the span of classical symbols (Li₅ at 13 S-units, all product shuffles incl. Li₃·Li₂ and Li₂·Li₂·log, log⁵) plus \(\mathrm{Li}_{3,2}\)-atom symbols (rank 32/32, consistent). One exact solution has 32 nonzero terms (6 Li₅, 9 Li₄·log, 7 Li₃·log², 4 Li₂·log³, 2 Li₃·Li₂, 2 Li₂²·log, coefficient height ~ \(3\times10^7\)) and the depth-two correction
\[
C \;=\; (-3600+2160\sqrt2)\,\mathrm{Li}_{3,2}(b,1)\;+\;(11880-8280\sqrt2)\,\mathrm{Li}_{3,2}(1/\sqrt2,1),
\]
whose pair equals the obstruction **exactly**; the classical remainder pairs to \((0,0)\). So \(S(W_5^{oct})=S(\text{classical}+C)\) is an exact symbol identity, and \(C\) is the octagon's own depth-two correction at M4c's arguments (not ratio·\(N_5\)).
- **NUMERICAL, inconclusive so far:** the residual \(R=W_5^{oct}-(\text{classical}+C)=2948.9558\ldots\) is exactly symbol-invisible, but LLL residual fits at dps 250 and dps 430 (detection threshold \(\sim10^{15}\)) over the 13-atom invisible basis \(\{\zeta5,\zeta2\zeta3,\zeta3\log^2(3),\zeta2\log^3(4),\zeta2^2\log(2),\zeta(3,2),K\}\) (\(K\) = the 1-dim symbol-kernel Li₃,₂ combination \(-8(b,1)-32(-b,1)-16(b,-1)-16(-b,-1)+(-b^2,1)\)) returned only LLL-threshold candidates — **no accepted relation**. The RREF solution is arbitrary within a 613-dimensional exact column-kernel, which can inflate residual coefficients well past \(10^{15}\) (M4c's analogous \(R_5\) needed 1200–1500 digits and had height \(1.4\times10^{15}\)); a kernel-combination completeness test was in progress when computation stopped. The classical-form hunt is therefore **not closed**; no claim is made either way beyond the PROVED symbol identity.

## Labels

- **PROVED** (exact symbol algebra, no floats): the rewrite algebra, all controls, the obstruction pairs, non-proportionality, the symbol-level identity \(S(W_5^{oct})=S(\text{classical}+C)\), \(C\)'s pair = obstruction.
- **NUMERICAL** (digits): all value-level checks listed in the evidence table, the \(G\)-table, and any closed-form statement (residual fit inconclusive).

## Evidence table

| # | Check | Result | Residual / agreement | Digits | Producer |
|---|---|---|---:|---:|---|
| 1 | 110-digit block targets from parsed_mons+itab110 | reproduced | sum−total \(7.7\times10^{-120}\); mixed−quoted \(-2.6\times10^{-50}\) | 120 | step_01.py |
| 2 | Exact rewrite of 756 integrals (1024 words, wt ≤ 5) | done; no leading-1 words | exact | — | oct5d2_engine.py, step (rewrite_exact.pkl) |
| 3 | G-table validation (polylogs, quad, nested quad, shuffles, reg values) | PASS | \(10^{-147}\)–\(10^{-171}\) (dps 170); \(10^{-252}\) (dps 250); \(10^{-432}\) (dps 430) | 170/250/430 | gengine.py |
| 4 | Rewrite: mixed block vs target/quoted | PASS | \(2.6\times10^{-107}\) / \(-2.6\times10^{-50}\) | 107/49 | step files, verify_oct5d2.py |
| 5 | Rewrite: single + local blocks | PASS | \(3.3\times10^{-106}\), \(2.8\times10^{-102}\) | 106/102 | same |
| 6 | Dual rewrite (\(u=1-t\)) numeric | PASS | \(1.2\times10^{-80}\) worst (40 random) | 80 | step files |
| 7 | Scaling invariance of symbols, 780 convergent words | PASS | exact, 0 mismatches | — | step files |
| 8 | Controls: Li₅(7 S-units), 300 shuffle products | pair to 0 | exact | — | verify_oct5d2.py |
| 9 | M4c W5 pairs to (84,−38√2), two implementations | PASS | exact | — | verify_oct5d2.py |
| 10 | G-symbol identities on test words | PASS | exact | — | verify_oct5d2.py |
| 11 | N₅ atom combination pairs to (84,−38√2) | PASS | exact | — | step files |
| 12 | Direct vs dual top symbols, all 756 integrals | identical | exact | — | step files |
| 13 | **Pairings: local 0; single-p = single-q; mixed; total \((-1485+1035\sqrt2, 900-540\sqrt2)\)** | **nonzero ⇒ depth two** | exact; dual assembly identical | — | verify_oct5d2.py |
| 14 | Non-proportionality to M₄ | ratios \(15(-33+23\sqrt2)/28 \ne 45(6-5\sqrt2)/19\) | exact | — | verify_oct5d2.py |
| 15 | Symbol span: \(S(W_5)\) = classical + Li₃,₂ symbols | rank 32/32, consistent; explicit 32-term solution | exact | — | symbol_fit2*.pkl, step files |
| 16 | Correction \(C\) pair = obstruction; classical part pair = 0 | PASS | exact | — | step files |
| 17 | Real-branch Liₙ(p) = Re(polylog), n=2..5 | PASS | \(0.0\) vs depth-1 integrals | 60 | step files |
| 18 | Residual invisible fit (dps 250, 430) | threshold only, not accepted | best maxc ~ \(6\times10^{14}\) at threshold ~ \(10^{15}\) | 430 | step files |

## Files in the work dir `new_results/OCT5D2_octagon_depth_test`

- **`verify_oct5d2.py`** — deliverable verifier: rewrite check (≥ 50 digits; achieves 83–86), G-table self-checks, all controls, all pairings incl. total and non-proportionality. Runs in **9.8 s**, prints 21 PASS lines, exit 0. (Verified with interpreter `<project>/.venv/bin/python`.)
- `step_01.py` — 110-digit block targets from the read-only inputs.
- `oct5d2_engine.py` — exact Q(√2) shuffle + IBP rewrite engine (Task 1 algebra).
- `gengine.py` — numerical G(w;1) tables (series at 1/2 + Chen; caches `gtab_d{dps}_N{N}.pkl`).
- `symengine.py` — symbol wrapper around the validated `A3_W4_proof/w4_symbol.py`.
- `result.md` — status report (includes all headline results).
- Data: `rewrite_exact.pkl` (exact rewrites), `block_symbols.pkl` (exact block+total weight-5 symbols), `targets110.pkl`, `w5_250.pkl`/`w5_430.pkl` (W₅/u to 240/420 digits), `atom_vals*.pkl`, `symbol_fit2*.pkl`, `K_val.pkl`, `kernel_sym.pkl` (exact 613-dim symbol-kernel basis for the ongoing invisible-completeness test).

**Provenance of reused validated code:** `qsqrt2.py` (paper `squircle_note/scripts`), `w4_symbol.py` (`A3_W4_proof`), M4c's `symbol.json`/`cobracket.json` — used exactly as in the M4c verification, ensuring the pairing functional is calibrated against M₄'s \((84,-38\sqrt2)\).
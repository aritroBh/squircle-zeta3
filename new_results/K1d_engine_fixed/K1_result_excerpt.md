The full 40-entry J-table (`Jtab_t.pkl`) was built symbolically by `esum_engine.py` (exact partial-fraction reduction of the Stirling/harmonic series $J=c!\sum_k s(k,c)\kappa^k M(n+k;a,b)/k!$) and **validated against direct quadrature at 40 tuples × 5 κ-values: worst residual 1.0e-31** (`step_31.py`).

### 2. The corrected general law (form derived; final numeric assembly not run)

For an edge with local term values $A_0t,\ A_1(1-t)$, one non-local term $B(1+\kappa t)$, $\kappa=r-1$ (so $A_0=A_1-B\kappa$, $A_1+B=1$), the $d^3$ integrand $E_3$ is the explicit $\mathbb Q$-polynomial expression of `E3_general.pkl` (40 non-local monomials $t^nX^aY^bW^c$ with coefficient functions of $(A_0,A_1,B,\kappa,\log A_i,\log B,L,q_0,|V_0|^2)$). Then
$$\boxed{\;P'''(1)\big|_{\zeta(3)}=\sum_e\bigl(Z_e+N_e(r)\bigr),\qquad N_e(r)=6\!\!\sum_{(n,a,b,c),\,c\ge1}\!\!\!\mathrm{coef}_{nabc}\;J(n,a,b,c;r-1)\;}$$
with $Z_e=6\alpha\beta(h^2/L)(3(\alpha u+\beta v)-2)$ (proved in E1) and $J$ in closed form via §1. For several non-local ratios per edge, cross terms $W_jW_k$ appear and need two-variable atoms (not yet built).

### 3. Special values making the counterexample decidable
- $A(-1)=-\tfrac58\zeta(3)$ (80-digit agreement; classical Euler) — governs $r=2$;
- $A(1/2)=\zeta(3)-\tfrac{\pi^2}{12}\log 2$ (PSLQ, residual $8.4\mathrm{e}{-81}$) — governs $r=\tfrac12$;
- $\mathrm{Li}_3(-1)=-\tfrac34\zeta(3)$, $\mathrm{Li}_3(1/2)=\tfrac78\zeta(3)-\tfrac{\pi^2}{12}\log2+\tfrac16\log^32$.

**HONEST FAILURE:** the final substitution of the edge coefficients into the J-table, which must output $N_e(2)=159\sqrt5/500$ per $r=2$ arc, plus the sheared case, the fresh 100-digit example, and the octagon decision, were **not executed** — the turn budget ended right after `E3_general.pkl` was produced. The check is fully assembled mechanically (coefficients × `Jtab_t.pkl` at $x=-1$ and $x=1/2$ with the values above) but I will not claim the number without printing it.

### 4. Octagon (partial, hand-derived + solver data)
Per edge of the regular octagon the two non-local ratios are $r=\sqrt2$ and $r=1/\sqrt2$ (computed from the vertices $(0,\sqrt2-1)$, $((1-\sqrt2)/\sqrt2,(\sqrt2-1)/\sqrt2)$). The atoms then sit at $x=1-\sqrt2$ and $x=1-\frac1{\sqrt2}=\frac{1-\sqrt2}{-\sqrt2}\cdot$…, i.e. related by the Landen map $x\mapsto x/(x-1)$, so three-term polylog identities may collapse the total. The solver's non-local total is NL $=0.08561735999445921640923680788284790860302$ (40 digits, `s13_out_40.json`); whether it equals a rational multiple of $\zeta(3)$ plus elementary terms — i.e. whether the octagon's $\zeta(3)$ coefficient is $C_4$ or $C_4+q\,\zeta(3)$ — is **still OPEN**: I did not run the PSLQ.

## Evidence table

| claim | status | evidence (file, residual) |
|---|---|---|
| my E3 reproduces verifier t5 decomposition | NUMERICAL | `step_12.py`, residual 5.0e-50 |

wd = '<project>/research/08_theory/K1f_hexagon_order5'
log_add = """step_24 (recovered run): rebuilt general single-ratio E5 TPS engine (approach A: Fpow=F^{-1/(1+d)}, metric L2-2yyp*Phi+yy*Phi^2, Phi=F'/F no 1/(1+d)); E[k] split into (t,X,Y,W) monomial dicts, saved E5_polys.pkl. E[3] vs K1d E3_fixed: 80 monomials, 0 mismatches (exact). E[5]: 336 monomials t^n X^a Y^b W^c, n<=5, a+b+c<=5, c>=1: 35 (a,b,c) patterns x n=0..5 = 210 nonlocal J-keys; local: 21 (a,b) x n=0..5 = 126 keys.
step_17/18 note: step_17 (E3 match) re-ran OK; step_18 (symb eval at t=0.3) TIMED OUT - replaced by plan: fast eval from E5_polys + jet check (to redo).
PLAN: (1) local Beta table Jloc(N;a,b) a+b<=5 exact in {1,z2,z3,z2^2,z2z3,z5} (E5 step_03 machinery, H_{N+1},H_{N+2} symbols); (2) esum engine weight 5 with H^{(1..4)} numerators, atoms raw symbols -> MPL -> HPL Lyndon basis {Li2..Li5, H011,H0011,H0111, H00011,H00101,H00111,H01011,H01111} + products + zetas; (3) J(n;a,b,c)=c! Sum_{m>=c} (-1)^{m-c} kap^m e_{c-1}(m-1)/m * Jloc(n+m;a,b); (4) assemble Earc5=5!*sum co*J; (5) checks: kappa->0 local law (Z_e^5, Z_e^{2,3}), numerics b0/b1/fresh/irrational hexagons vs per-t Cauchy+tanh-sinh; (6) zeta5/zeta2zeta3 content.
"""
with open(f'{wd}/progress.log','a') as f:
    f.write(log_add)

result = r"""# K1f — P^(5)(1) of EVERY centrally symmetric hexagon in closed form (DRAFT, in progress)

Work dir: `research/08_theory/K1f_hexagon_order5/`. Weight-5 version of K1e. Reuses: K1d fixed E3 builder convention (validated), E5 local laws, K1c e-sum J-table approach.

## Conventions (locked, matching K1d/K1e numerics)
- Edge data: term values a0=A0 t, a1=A1(1-t), a2=B(1+kap t), kap=r-1, A1+B=1, A0=1-B r; L2=|V1-V0|^2, q0=V0.(V1-V0), V02=|V0|^2; unweighted logs lA0=log(A0/w0), lA1=log(A1/w1), lB=log(B/wo).
- X=log t, Y=log(1-t), W=log(1+kap t) (kept symbolic; the ONLY non-local log).
- Atom argument x = -kap = 1-r. E_arc^(5) = 5! * int_0^1 E5(t) dt; P^(5)(1) = 2 * Sum_{q=0}^2 E_arc^(5)(edge_q).
- L1 atom convention: L1 = -ln(1-x) = -log r (K1e planner correction respected).

## Status
- Task 1: general E5(t) BUILT exactly (step_24: TPS engine to d^5; E[3] matches K1d E3_fixed EXACTLY, 80 monomials 0 mismatches; E[5] = 336 monomials in (t,X,Y,W), n<=5). Numeric t=0.3 jet check: PENDING (fast eval path).
- Task 2: weight-5 J-table: IN PROGRESS (local table + e-sum extension + HPL Lyndon reduction).
- Task 3: assembly + checks: PENDING.
- Task 4: zeta(5), zeta(2)zeta(3) content: PENDING.
"""
with open(f'{wd}/result.md','w') as f:
    f.write(result)
print("result.md + progress.log written")

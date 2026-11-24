
import os
WD = '<project>/research/08_theory/K1f_hexagon_order5'
os.makedirs(WD, exist_ok=True)
with open(f'{WD}/progress.log','a') as f:
    f.write("step_01/02: explored research tree; read K1e result.md (incl planner correction on l1=-ln(1-x) sign), K1d step_33_builder (E3 general single-ratio builder), K1c esum_engine (harmonic-sum J-reduction), E1 s01_edge_sym (authoritative convention |x_p'|=F^{-1/p}sqrt((L-sig Phi)^2+h^2 Phi^2), P^{(k)}(1)=k! Sum_e int E_k), E5 step_01/02 (local TPS engine, tau/gam lumped), E5 result.md. Resolved F-exponent: K1e G=F^{-2/(1+d)}*metric is pre-sqrt; sqrt(G)=F^{-1/(1+d)}*sqrt(metric) = E1/E5 convention. Plan: build general single-ratio E5 by un-lumping tau/gam in E5 step_02 engine.\n")
print("progress.log updated")

# step_03: derive d^n coefficient of sqrt(sum g_i d^i), confirm d^3 matches K1d's corrected formula
import sympy as sp
d = sp.symbols('d')
g = sp.symbols('g0:6')
Ls = sp.sqrt(g[0])
# series of sqrt
ser = sp.sqrt(sum(g[i]*d**i for i in range(6)) + sp.O(d**6)).series(d, 0, 6).removeO().expand()
for k in range(6):
    ck = sp.expand(ser.coeff(d, k))
    print(f"[d^{k}] sqrt = {ck}")
# confirm d^3 coefficient equals g3/(2Ls) - g1*g2/(4*Ls^3) + g1^3/(16*Ls^5)
c3_known = g[3]/(2*Ls) - g[1]*g[2]/(4*Ls**3) + g[1]**3/(16*Ls**5)
print("d^3 matches K1d corrected formula:", sp.simplify(ser.coeff(d,3) - c3_known)==0)

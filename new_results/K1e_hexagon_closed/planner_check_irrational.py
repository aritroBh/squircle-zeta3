# Planner check: K1e formula (Earc_basis_clean.pkl, implementation convention l1 = -log r) vs direct polar computation,
# golden hexagon V = (1,0), (phi,1), (0,1) (+negatives), ratios phi, 1, 1/phi (irrational), weights w = (2,3,5).
import pickle, sympy as sp, mpmath as mp
mp.mp.dps = 60
WORK = 'research/08_theory/K1e_hexagon_closed'
phi = (1+sp.sqrt(5))/2
V = [sp.Matrix([1, 0]), sp.Matrix([phi, 1]), sp.Matrix([0, 1])]
w = [2, 3, 5]
cross = lambda P, Q: P[0]*Q[1]-P[1]*Q[0]
c = sp.symbols('c0:3')   # c_j = w_j k_j, from N(V_q) = 1
sol = sp.solve([sum(c[j]*sp.Abs(cross(V[j], V[q])) for j in range(3)) - 1 for q in range(3)], c, dict=True)[0]
k = [sp.nsimplify(sp.simplify(sol[c[j]]/w[j])) for j in range(3)]
print("k =", k)
clean = pickle.load(open(f"{WORK}/Earc_basis_clean.pkl", "rb"))
L2s = sp.Symbol('L2'); pars = sp.symbols('A1 B kap q0 V02 lA0 lA1 lB'); atoms = sp.symbols('l3a l3l l2l1 l1_3 l2a l1_2 l1 z3 z2l1 z2 one')
order = ['Li3(x)','Li3(x/(x-1))','Li2(x)*ln(1-x)','ln^3(1-x)','Li2(x)','ln^2(1-x)','ln(1-x)','zeta3','zeta2*ln(1-x)','zeta2','1']
f = sp.lambdify(pars + (L2s,) + atoms, sum(num/L2s**sp.Rational(5, 2)*a for (num, den), a in zip([clean[n] for n in order], atoms)), 'mpmath')
verts = [V[0], V[1], V[2], -V[0]]
tot = 0
for e in range(3):
    V0, V1 = verts[e], verts[e+1]
    j0 = [j for j in range(3) if sp.simplify(cross(V[j], V0)) == 0][0]
    j1 = [j for j in range(3) if sp.simplify(cross(V[j], V1)) == 0][0]
    jo = 3 - j0 - j1
    A0 = sp.simplify(w[j0]*sp.Abs(k[j0]*cross(V[j0], V1))); A1 = sp.simplify(w[j1]*sp.Abs(k[j1]*cross(V[j1], V0)))
    B = sp.simplify(w[jo]*sp.Abs(k[jo]*cross(V[jo], V0))); r = sp.simplify(cross(V[jo], V1)/cross(V[jo], V0))
    assert sp.simplify(A1 + B - 1) == 0 and sp.simplify(A0 - (1 - B*r)) == 0
    L2 = sp.simplify((V1-V0).dot(V1-V0)); q0 = sp.simplify(V0.dot(V1-V0)); V02 = sp.simplify(V0.dot(V0))
    x = 1 - mp.mpf(sp.N(r, 80)); l1 = -mp.log(1-x)
    val = f(mp.mpf(sp.N(A1, 80)), mp.mpf(sp.N(B, 80)), mp.mpf(sp.N(r-1, 80)), mp.mpf(sp.N(q0, 80)), mp.mpf(sp.N(V02, 80)),
            mp.log(mp.mpf(sp.N(A0/w[j0], 80))), mp.log(mp.mpf(sp.N(A1/w[j1], 80))), mp.log(mp.mpf(sp.N(B/w[jo], 80))), mp.mpf(sp.N(L2, 80)),
            mp.polylog(3, x), mp.polylog(3, x/(x-1)), mp.polylog(2, x)*l1, l1**3, mp.polylog(2, x), l1**2, l1, mp.zeta(3), mp.pi**2/6*l1, mp.pi**2/6, mp.mpf(1))
    print(f"edge {e+1}: r = {sp.nsimplify(r)} ~ {mp.nstr(mp.mpf(sp.N(r,30)),12)}, E_arc = {mp.nstr(val, 30)}")
    tot += val
formula = 2*tot
# direct polar computation
kk = [mp.mpf(sp.N(kj, 80)) for kj in k]; Vn = [(mp.mpf(sp.N(v[0], 80)), mp.mpf(sp.N(v[1], 80))) for v in V]
def P(p):
    cc = [w[j]*kk[j]**p for j in range(3)]
    def g(th):
        s, co = mp.sin(th), mp.cos(th); F = 0; Fp = 0
        for j in range(3):
            sj = Vn[j][0]*s - Vn[j][1]*co; dsj = Vn[j][0]*co + Vn[j][1]*s
            F += cc[j]*abs(sj)**p; Fp += cc[j]*p*abs(sj)**(p-1)*mp.sign(sj)*dsj
        rho = F**(-1/p); drho = -(1/p)*F**(-1/p-1)*Fp
        return mp.sqrt(rho**2 + drho**2)
    br = [mp.atan2(v[1], v[0]) for v in Vn] + [mp.pi]
    return 2*sum(mp.quad(g, [br[i], br[i+1]]) for i in range(3))
d3 = lambda h: (P(1+2*h) - 2*P(1+h) + 2*P(1-h) - P(1-2*h))/(2*h**3)
h = mp.mpf(10)**-6; a, b = d3(h), d3(h/2); direct = (4*b-a)/3
print("P(1) direct =", mp.nstr(P(mp.mpf(1)), 20))
print("formula P'''(1) =", mp.nstr(formula, 30)); print("direct  P'''(1) =", mp.nstr(direct, 30)); print("diff =", mp.nstr(formula-direct, 3))

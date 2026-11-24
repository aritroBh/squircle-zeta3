"""Numerical check (referee quadrature) of every J-table entry and every local integral used by K1e,
including kappa >= 1 where the kappa power series used to build the table diverges."""
import sympy as sp, mpmath as mp, pickle, json
mp.mp.dps = 45
D = '<project>/research/08_theory/K1e_hexagon_closed/'
J = pickle.load(open(D+'Jtab_t.pkl', 'rb'))
jl = json.load(open(D+'jloc.json'))
L1s, L2s, L3s, AAs, z2s, xs = sp.symbols('L1 L2 L3 AA z2 x')
def Acont(x):
    l = mp.log(1-x)
    return (12*mp.polylog(3,x) + 6*mp.polylog(3, x/(x-1)) - 6*mp.polylog(2,x)*l - l**3)/6
worst = {}
for kap in ['-0.9', '-0.5', '0.4', '1', '1.7', '3.5', '9']:
    k = mp.mpf(kap); x = -k
    vals = {L1s: -mp.log(1-x), L2s: mp.polylog(2,x), L3s: mp.polylog(3,x), AAs: Acont(x), z2s: mp.zeta(2), xs: x}
    w = 0
    for (n,a,b,c), e in J.items():
        f = sp.lambdify(list(vals), e, 'mpmath')
        cf = f(*vals.values())
        num = mp.quad(lambda t: t**n*mp.log(t)**a*mp.log(1-t)**b*mp.log(1+k*t)**c, [0, mp.mpf(1)/2, 1])
        w = max(w, abs(cf-num)/max(1, abs(num)))
    worst[kap] = w
    print('kappa', kap, ' worst rel. error over', len(J), 'J entries:', mp.nstr(w, 3), flush=True)
z2, z3 = sp.symbols('z2 z3'); w = 0
for key, v in jl.items():
    n, a, b = map(int, key.split(','))
    cf = sp.sympify(v, locals={'z2': z2, 'z3': z3}).subs({z2: sp.zeta(2), z3: sp.zeta(3)})
    num = mp.quad(lambda t: t**n*mp.log(t)**a*mp.log(1-t)**b, [0, mp.mpf(1)/2, 1])
    w = max(w, abs(mp.mpf(sp.N(cf, 60)) - num))
print('local Beta integrals (', len(jl), '): worst abs error', mp.nstr(w, 3))

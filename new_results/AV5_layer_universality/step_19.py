import sympy as sp, mpmath, pickle, json
from mpmath import mp, mpf, sqrt, zeta, log, exp, quad, tanh, nstr, factorial, polylog
mp.dps = 90

# ---- symbolic f_k with placeholder symbols g1,g3,g5,g7 ----
tau, s = sp.symbols('tau s')
G1,G3,G5,G7 = sp.symbols('G1 G3 G5 G7')
Gser = tau*G1 + tau**3*G3 + tau**5*G5 + tau**7*G7
eG = sp.series(sp.exp(Gser), tau, 0, 9).removeO().expand()
css = sp.series(sp.cos(tau*s)+sp.sin(tau*s), tau, 0, 9).removeO().expand()
prod = sp.expand((eG-1)*css)
Fk = {k: sp.expand(prod.coeff(tau,k)) for k in range(1,9)}
# each Fk is a polynomial in s, G1,G3,G5,G7. Extract monomials s^a G1^b G3^c G5^d G7^e
mons = {k: [] for k in range(1,9)}
for k,expr in Fk.items():
    for term in sp.Add.make_args(expr):
        coeff, rest = term.as_coeff_Mul()
        powers = {G1:0,G3:0,G5:0,G7:0,s:0}
        for fct in sp.Mul.make_args(rest):
            base, ex = fct.as_base_exp()
            powers[base] = ex
        mons[k].append((sp.Rational(coeff), powers[s], powers[G1], powers[G3], powers[G5], powers[G7]))
for k in range(1,9):
    print(k, len(mons[k]), "monomials; e.g.", mons[k][:3])
pickle.dump({k:[(str(c),a,b,cc,d,e) for (c,a,b,cc,d,e) in v] for k,v in mons.items()},
            open('ckpt_monomials.pkl','wb'))

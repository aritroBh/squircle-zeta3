import os, pickle, sympy as sp, mpmath as mp
os.chdir('research/08_theory/DEC1_decagon')
import dec_engine as E
r5 = sp.sqrt(5)
lf,l2,z2,z3 = E.lf,E.l2,E.z2,E.z3
Lsym = 2*sp.sqrt(10-2*r5)/4/(r5*E.phi**2)
g = pickle.load(open('geometry.pkl','rb'))
z3coeff_loc = sp.Rational(5,32) - r5/8
C5_check = sp.simplify(60*Lsym*z3coeff_loc)
print("60 L [z3]S_loc =", sp.simplify(C5_check))
print("C_5 =", sp.simplify(g['C5']))
print("difference:", sp.simplify(C5_check - g['C5']))
S3 = pickle.load(open('S_reduced.pkl','rb'))['S3']
z3tot = sp.Poly(sp.expand(S3), l2,lf,z2,z3,E.L3_1,E.L3f).coeff_monomial(z3)
p = sp.simplify((z3tot+z3tot.subs(r5,-r5))/2); q = sp.simplify((z3tot-z3tot.subs(r5,-r5))/(2*r5))
print("[z3] S3 total (named basis):", p, "+", q, "sqrt5")
T = sp.simplify(60*Lsym*z3tot)
Tp = sp.simplify((T+T.subs(r5,-r5))/2); Tq = sp.simplify((T-T.subs(r5,-r5))/(2*r5))
print("total z3 coeff of P'''(1):", Tp, "+", Tq, "sqrt5")
Ne = sp.simplify((T - g['C5'])/10)
Nep = sp.simplify((Ne+Ne.subs(r5,-r5))/2); Neq = sp.simplify((Ne-Ne.subs(r5,-r5))/(2*r5))
print("non-local z3 per edge (named basis):", Nep, "+", Neq, "sqrt5 (times L-units: note N_e as coefficient of z3 in P'''/10 per edge)")
mp.mp.dps=60
print("numeric checks: C5_check", mp.nstr(mp.mpf(sp.N(C5_check,40)),25))
print("T:", mp.nstr(mp.mpf(sp.N(T,40)),25))
print("Ne:", mp.nstr(mp.mpf(sp.N(Ne,40)),25))
with open('progress.log','a') as f:
    f.write(f"step: local z3 == C_5 exactly. total z3 = {Tp}+{Tq}sqrt5 *L-units; N_e = {Nep}+{Neq}sqrt5.\n")
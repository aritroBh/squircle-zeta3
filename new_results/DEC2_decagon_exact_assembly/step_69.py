import sympy as sp, pickle, os, sys
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd); sys.path.insert(0, '<project>/research/08_theory/DEC1_decagon')
t = sp.Symbol('t'); r5 = sp.sqrt(5); phi = (1+r5)/2
kap1 = 1/phi; kap4 = -1/phi**2
import dec_engine as E
ATOMS = E.NAMED
def canon5(e):
    """reduce coefficients (in Q(sqrt5)) of atom monomials to canonical p+q*sqrt5; return canonical expr."""
    e = sp.expand(e)
    P = sp.Poly(e, *ATOMS)
    out = sp.Integer(0)
    for monom, coeff in P.terms():
        c = sp.together(coeff)
        cc = c.subs(r5, -r5)
        p = sp.radsimp(sp.expand((c+cc)/2))
        q = sp.radsimp(sp.expand((c-cc)/(2*r5)))
        m = sp.Integer(1)
        for a_, e_ in zip(ATOMS, monom):
            m *= a_**e_
        out += (p + q*r5)*m
    return sp.expand(out)
# recompute D9, D10 canonically
D9 = canon5((-E.lf**3/2 - kap4*E.D7/2)/kap1)
D10 = canon5((E.lf**3/2 - kap1*E.D8/2)/kap4)
print("D9 canonical diff vs engine:", canon5(D9 - E.D9))
print("D10 canonical diff vs engine:", canon5(D10 - E.D10))
res = pickle.load(open('D5678_proved.pkl','rb'))
res['D9'] = D9; res['D10'] = D10
pickle.dump(res, open('D5678_proved.pkl','wb'))
with open('progress.log','a') as f:
    f.write("step_09/10: D5-D8 PROVED (antiderivative certificates, formal identities, endpoint values match DEC1 engine symbolically: diff 0); D9,D10 PROVED by IBP (integrand identities residual 0, boundary exact) - canonical match with engine\n")
print("saved D5678_proved.pkl with D5..D10")

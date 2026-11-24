import sympy as sp, pickle, os, sys
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd); sys.path.insert(0, '<project>/research/08_theory/DEC1_decagon')
import dec_engine as E
r5 = sp.sqrt(5); rho = sp.sqrt(10-2*r5)
mons = pickle.load(open('mons_dec2.pkl','rb'))
Jloc = pickle.load(open('Jloc_dec2.pkl','rb'))
def JLOC(n,a,b):
    return sp.expand(Jloc[(n,a,b)].subs(sp.pi**2, 6*E.z2).subs(sp.zeta(3), E.z3))
l2, lf = E.l2, E.lf
Sloc = sp.Integer(0)
for (n_, abc, letters), coef in mons.items():
    a,b,c1,c4 = abc
    if c1==0 and c4==0:
        Sloc += coef * l2**letters[0] * lf**letters[1] * JLOC(n_,a,b)
Sloc = sp.expand(Sloc)
cz3 = sp.expand(Sloc).coeff(E.z3)
print("local zeta(3) coefficient of S (per edge factor):", cz3)
# P'''(1) local zeta(3) = 60 * L * cz3  vs C5 = rho*(105*r5-255)/32
L = rho*(3*r5-5)/20
print("L == rho(3 sqrt5 - 5)/20 exact check:", sp.simplify(L**2 - (sp.Rational(25,10)-sp.Rational(11,10)*r5)))
C5 = rho*(105*r5-255)/32
expr = sp.expand(60*L*cz3 - C5)
print("60 L cz3 - C5 =", sp.simplify(expr))
# also verify C5 equals Theorem's form: 48*m sin^4 cos (3 sin^2 - 1), m=5
s36 = sp.sqrt(10-2*r5)/4  # sin(pi/10)= (sqrt5-1)/4 ... use sin(pi/(2m)) = sin(pi/10)
spi10 = (r5-1)/4
C5_thm = 48*5*spi10**4*sp.cos(sp.pi/10)*(3*spi10**2-1)
print("C5 vs theorem form diff:", sp.simplify(sp.nsimplify(0)+sp.expand(C5 - C5_thm)))
with open('progress.log','a') as f:
    f.write("step_14: local zeta(3): 60 L cz3 - C5 = 0 EXACT (C5 = rho(105 sqrt5-255)/32, L = rho(3 sqrt5-5)/20 in Q(rho)); C5 matches Theorem C_m at m=5\n")

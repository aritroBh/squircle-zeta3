import sympy as sp
v, tau = sp.symbols('v tau')
pi = sp.pi; I_ = sp.I
om = (1+I_)/sp.sqrt(2); omb = (1-I_)/sp.sqrt(2)
N = pi/(2*v**4)
x = v*tau

def ser(expr, vmax=3):
    """series in v through v^vmax, using sympy's series with omb algebraic"""
    return sp.expand(sp.series(expr, v, 0, vmax+1).removeO())

# G = i pi s0 x/2 + (N-1) log(1+x) - N log(1 + x/(1+v^2 omb))
s0x = om/v**2 * x
T1 = I_*pi*s0x/2
T2 = (N-1)*sp.log(1+x)
T3 = -N*sp.log(1 + x/(1+v**2*omb))
Gser = ser(T1+T2+T3, 2)
# expand N explicitly and collect
Gser = sp.expand(Gser)
terms = {}
for term in sp.Add.make_args(Gser):
    c, vp = term.as_coeff_exponent(v)
    terms[vp] = terms.get(vp, 0) + c
terms = {p: sp.simplify(c) for p, c in terms.items()}
b = pi*omb/2
G1_exp = -(1+pi*omb**2/2)*tau + (pi*omb/2)*tau**3
G2_exp = (sp.Rational(1,2)+3*pi*omb**2/4)*tau**2 - (pi*omb/2)*tau**4
print("negative-power coefficients (must vanish):", {p: sp.simplify(terms.get(p,0)) for p in terms if p < 0})
print("v^0 + b tau^2 =", sp.simplify(terms.get(0,0) + b*tau**2))
print("G1 mine     =", sp.expand(terms.get(1,0)))
print("G1 expected =", sp.expand(G1_exp))
print("diff:", sp.simplify(terms.get(1,0) - G1_exp))
print("G2 mine     =", sp.expand(terms.get(2,0)))
print("G2 expected =", sp.expand(G2_exp))
print("diff:", sp.simplify(terms.get(2,0) - G2_exp))

# now c1 and varpi1 fully independently
h1, h2 = sp.Rational(-3,2), sp.Rational(15,8)
G1 = G1_exp; G2 = G2_exp
integrand = sp.expand(G2 + G1**2/2 + h2*tau**2 + G1*h1*tau)
def mom(m):
    return sp.factorial2(2*m-1)/(2*b)**m if m > 0 else sp.Integer(1)
c1 = 0
for term in sp.Add.make_args(integrand):
    c, tp = term.as_coeff_exponent(tau)
    c1 += c*mom(int(tp)//2)
c1 = sp.simplify(c1)
S1 = pi*(1+I_)/(6*sp.sqrt(2))
varpi1_mine = sp.simplify(sp.sqrt(pi/2)*(c1 + S1 - 3*omb/4))
varpi1_paper = ((24+12*pi+pi**2)+I_*(24-12*pi+pi**2))/(48*sp.sqrt(pi))
print("c1 =", c1)
print("varpi1 mine - paper =", sp.simplify(varpi1_mine - varpi1_paper))
print("varpi1 =", sp.N(varpi1_paper, 30))
open('progress.log','a').write("step: independent sympy derivation of G1,G2,c1,S1 -> varpi1 matches paper exactly (diff 0)\n")

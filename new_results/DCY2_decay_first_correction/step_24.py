import sympy as sp
pi = sp.pi; I_ = sp.I
wb = (1-I_)/sp.sqrt(2)
b, tau = sp.symbols('b tau')
bval = pi*wb/2

G1 = -(1 + pi*wb**2/2)*tau + (pi*wb/2)*tau**3
G2 = (sp.Rational(1,2) + 3*pi*wb**2/4)*tau**2 - (pi*wb/2)*tau**4
h1 = sp.Rational(-3,2); h2 = sp.Rational(15,8)

integrand_v2 = sp.expand(G2 + G1**2/2 + h2*tau**2 + G1*h1*tau)

def mom(m):  # <tau^{2m}> for exp(-b tau^2)
    return sp.ff(2*m-1, m) / (2*b)**m if m > 0 else 1

res = 0
for term in sp.Add.make_args(integrand_v2):
    c, tpow = term.as_coeff_exponent(tau)
    res += c*mom(int(tpow)//2)
c1 = sp.simplify(sp.simplify(res.subs(b, bval)).rewrite(sp.exp).expand())
print("c1 =", c1)

S1 = pi*(1+I_)/(6*sp.sqrt(2))
pref = S1 - 3*wb/4
varpi1_mine = sp.simplify(sp.sqrt(pi/2)*(c1 + pref))
varpi1_paper = ((24+12*pi+pi**2) + I_*(24-12*pi+pi**2))/(48*sp.sqrt(pi))
print("mine  =", sp.N(varpi1_mine, 30))
print("paper =", sp.N(varpi1_paper, 30))
print("diff  =", sp.simplify(varpi1_mine - varpi1_paper))

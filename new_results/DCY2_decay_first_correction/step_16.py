import sympy as sp

# Exact symbols: pi, and wb = exp(-i pi/4) = (1-i)/sqrt(2)
pi = sp.pi
I_ = sp.I
wb = (1-I_)/sp.sqrt(2)
b, tau = sp.symbols('b tau')

# b = pi*wb/2
bval = pi*wb/2

# G1, G2 coefficients (functions of tau)
G1 = -(pi*wb**2/2)*tau + (pi*wb/2)*tau**3
G2 = (3*pi*wb**2/4)*tau**2 - (pi*wb/2)*tau**4
h1 = sp.Rational(-3,2); h2 = sp.Rational(15,8)

integrand_v2 = G2 + G1**2/2 + h2*tau**2 + G1*h1*tau
integrand_v2 = sp.expand(integrand_v2)

# Gaussian moments <tau^{2m}> wrt exp(-b tau^2): (2m-1)!!/(2b)^m
def mom(m):
    return sp.ff(2*m-1, m) / (2*b)**m if m>0 else 1

# expectation
expr = integrand_v2
res = 0
for term in sp.Add.make_args(expr):
    c, tpow = term.as_coeff_exponent(tau)
    m = int(tpow)//2
    res += c*mom(m)
res = sp.simplify(res)
# substitute b
c1 = sp.simplify(res.subs(b, bval))
c1 = sp.simplify(c1.rewrite(sp.exp))
print("c1 =", sp.simplify(c1))

# varpi1 = sqrt(pi/2) * c1 ; compare with paper's closed form
varpi1_paper = ((24+12*pi+pi**2) + I_*(24-12*pi+pi**2))/(48*sp.sqrt(pi))
varpi1_mine = sp.simplify(sp.sqrt(pi/2)*c1)
print("mine =", sp.simplify(varpi1_mine))
print("diff =", sp.simplify(varpi1_mine - varpi1_paper))
print("numeric mine:", sp.N(varpi1_mine, 30))
print("numeric paper:", sp.N(varpi1_paper, 30))

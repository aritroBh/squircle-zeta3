import sympy as sp

v, tau = sp.symbols('v tau')
pi = sp.pi; I_ = sp.I
om = (1+I_)/sp.sqrt(2); omb = (1-I_)/sp.sqrt(2)
N = pi/(2*v**4)
s0 = om/v**2

# phi(s) = i pi s/2 + (N-1) log s - N log(1+s);  s = s0 (1 + v tau)
# G = phi(s0(1+v tau)) - phi(s0)
w = sp.symbols('w')  # local variable for series
x = v*tau
# term1: i pi s0 x /2
T1 = I_*pi*s0*x/2
# term2: (N-1) log(1+x)
ser_log1 = sum((-1)**(j+1)*x**j/j for j in range(1,6))
T2 = (N-1)*ser_log1
# term3: -N log(1 + s0 x/(1+s0)),  s0/(1+s0) = 1/(1+y omb), y = v^2
q = x/(1+v**2*omb)
T3 = -N*sum((-1)**(j+1)*q**j/j for j in range(1,6))
G = sp.expand(T1+T2+T3)
# collect by powers of v
Gser = sp.Poly(sp.expand(G), v)
terms = {}
for (p,), c in zip(Gser.terms(), Gser.coeffs()):
    terms[p] = sp.simplify(c)
for p in sorted(terms):
    if p <= 2:
        print(f"v^{p}: ", sp.simplify(terms[p]))
# expected
b = pi*omb/2
G1_exp = -(1+pi*omb**2/2)*tau + (pi*omb/2)*tau**3
G2_exp = (sp.Rational(1,2)+3*pi*omb**2/4)*tau**2 - (pi*omb/2)*tau**4
print("const term (should be -b tau^2):", sp.simplify(terms.get(0,0) + b*tau**2))
print("G1 diff:", sp.simplify(terms.get(1,0) - G1_exp))
print("G2 diff:", sp.simplify(terms.get(2,0) - G2_exp))

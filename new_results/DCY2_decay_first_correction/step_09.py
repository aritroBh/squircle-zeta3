import sympy as sp

# Symbols
pi = sp.pi
N = sp.symbols('N', positive=True)
w = sp.Symbol('w')  # perturbation s = s0 + w

# s0 = rho*omega as exact algebraic: rho = sqrt(2)*N**(1/2)/sqrt(pi); omega = (1+i)/sqrt(2)
i = sp.I
s0 = sp.sqrt(N)/(sp.sqrt(pi)) * (1+i)   # rho*omega = sqrt(2N/pi)*(1+i)/sqrt2 = sqrt(N/pi)(1+i)

# Exact phi'(s) = (N-1)/s - N/(1+s) + i*pi/2 ; series at s0 through w^4
s = sp.Symbol('s')
phi1 = (N-1)/s - N/(1+s) + i*pi/2
ser = sp.series(phi1, s, s0, 5).removeO().expand()
# Extract coefficients of w = s - s0
poly = sp.Poly(sp.expand(ser), w) if False else None
# series gives in (s - s0); replace
expr = sp.expand(ser)
# Use t = s - s0
t = sp.Symbol('t')
expr_t = sp.expand(expr.subs(s, t + s0))
c = [sp.simplify(expr_t.coeff(t, j)) for j in range(5)]
for j, cj in enumerate(c):
    print(j, sp.simplify(cj))

# SM5 method test: squircle c5 from Theorem 1 at N=5 (should be -sqrt2(3/4 z5 + 7/4 z2 z3)).
import sympy as sp
s, t, tau = sp.symbols('s t tau', positive=True)
Y = sp.Symbol('Y'); E = sp.Symbol('E')
# Lambda(x) = artanh(tan x): series, then y = 2 Lambda(tau s)/tau
x = sp.Symbol('x')
Lam = sp.series(sp.atanh(sp.tan(x)), x, 0, 8).removeO()
y = sp.expand(2*Lam.subs(x, tau*s)/tau)
y = sum(y.coeff(tau, k)*tau**k for k in range(0, 5))
print("y =", y)
d = sp.expand(y - 2*s)
Fgen = sp.log(1+sp.exp(-Y))
Ls = sp.expand(tau*sum(sp.diff(Fgen, Y, k).subs(Y, 2*s)*d**k/sp.factorial(k) for k in range(4)))
g = {k: sp.simplify(Ls.coeff(tau, k)) for k in range(0, 5)}
print({k: g[k] for k in range(1, 5)})
assert g[2] == 0 and g[4] == 0
g1, g3 = g[1], g[3]
gp = sp.diff(g1, s)
kappa = sp.simplify(g3/(s**3*gp)); print("g3 = kappa s^3 g1' with kappa =", kappa)
H = 1/sp.sqrt(2)
gg1, gg3 = sp.symbols('g1 g3')
gen = sp.series(sp.exp(t*gg1 + t**3*gg3) - 1, t, 0, 5).removeO()
G = {n: sp.expand(gen.coeff(t, n)) for n in range(1, 5)}
Mgen = sp.series(H*sp.cos(t*s) + H*sp.sin(t*s), t, 0, 4).removeO()
M = {a: sp.expand(Mgen.coeff(t, a)) for a in range(4)}
Phi4 = sp.expand(M[0]*G[4] + M[1]*G[3] + M[2]*G[2] + M[3]*G[1])
GP = sp.Symbol('gp'); Gs = sp.Symbol('g')
Phi4 = sp.expand(Phi4.subs({gg3: kappa*s**3*GP}).subs(gg1, Gs))
gm = lambda a, b: sp.Symbol('gam_%d%d' % (a, b))
total = 0
for term in sp.Add.make_args(Phi4):
    c, rest = term.as_independent(s, Gs, GP)
    pw = rest.as_powers_dict(); k, j, e = pw.get(s, 0), pw.get(Gs, 0), pw.get(GP, 0)
    if e == 1: c, k, j = -c*sp.Rational(k, j+1), k-1, j+1
    total += c*sp.factorial(k)*sp.factorial(j)*gm(k, j)
total = sp.expand(total); print("int_0^inf Phi4 =", total)
z2, z3, z5 = sp.symbols('z2 z3 z5')
rel = {gm(0,4): -gm(1,3) + 2*gm(3,1) - z2*z3/8 + z5/8,      # from sm5_derive.py (Lemma 5, degree 4)
       gm(2,2): -3*gm(3,1) + z2*z3/16 + z5/16}
val = sp.expand(total.subs(rel))
g31 = sp.Rational(15, 256)*z5                             # gamma_31 = (1/6)(3/8)eta(5)
print("after reflection relation:", val)
res = sp.expand(8*val.subs(gm(3,1), g31))
print("c5(squircle) =", res)
target = -sp.sqrt(2)*(sp.Rational(3,4)*z5 + sp.Rational(7,4)*z2*z3)
print("difference from -sqrt2(3/4 z5+7/4 z2 z3):", sp.simplify(res - target))

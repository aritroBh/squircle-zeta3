# SM5 item 2: exact derivation of c5 for the soft-max from Theorem 1 at N=5 (sympy).
import sympy as sp
s, t = sp.symbols('s t', positive=True)
E = sp.Symbol('E', positive=True)          # E = e^{-2s}
# --- step 1: layer expansion.  L = tau*F(y), y = 2 sin(tau s)/tau, F(y)=log(1+e^{-y})
tau = sp.Symbol('tau', positive=True)
y = sp.series(2*sp.sin(tau*s)/tau, tau, 0, 6).removeO()
print("y =", sp.expand(y))
d = sp.expand(y - 2*s)                      # shift delta
Fy = sp.log(1+E)                            # F(2s)
F1 = -1/(1+1/E)                             # F'(y)   = -1/(1+e^y), e^y = 1/E at y=2s
F2 = sp.diff(F1, E)*(-2*E)                  # d/dy = (dE/dy) d/dE, dE/dy = -E  (E=e^{-y}); here y-derivative -> use E=e^{-y}
# do it cleanly: F as function of y via E=e^{-y}
Y = sp.Symbol('Y')
Fgen = sp.log(1+sp.exp(-Y))
Ls = tau*sum(sp.diff(Fgen, Y, k).subs(Y, 2*s)*d**k/sp.factorial(k) for k in range(0, 4))
Ls = sp.expand(Ls)
g = {k: sp.simplify(Ls.coeff(tau, k)) for k in range(0, 7)}
for k in range(1, 6): print("g%d ="%k, g[k])
assert g[2] == 0 and g[4] == 0
g1, g3, g5 = g[1], g[3], g[5]
# --- step 2: generating functions (A=1, so P_b=delta_{b0}); tent for s>=0: H cos + T sin, H=T=1
gg1, gg3 = sp.symbols('g1 g3')             # abstract profile values
gen = sp.series(sp.exp(t*gg1 + t**3*gg3) - 1, t, 0, 5).removeO()
G = {n: sp.expand(gen.coeff(t, n)) for n in range(1, 5)}
Mgen = sp.series(sp.cos(t*s) + sp.sin(t*s), t, 0, 4).removeO()
M = {a: sp.expand(Mgen.coeff(t, a)) for a in range(0, 4)}
print("M =", M); print("G =", G)
Phi4 = sp.expand(M[0]*G[4] + M[1]*G[3] + M[2]*G[2] + M[3]*G[1])
print("Phi4 =", Phi4)
# --- step 3: reduce by parts.  g3 = -(s^3/6) g1'  with g1' = -2/(1+e^{2s})
gp = sp.Symbol('gp')                        # stands for g1'
assert sp.simplify(g3 - (-(s**3/6)*sp.diff(g1, s))) == 0
Phi4_g = sp.expand(Phi4.subs(gg3, -(s**3/6)*gp))
print("Phi4 in (s, g, g') =", Phi4_g)
# monomials s^k g^j gp^e (e in {0,1}); integrate over (0,inf), boundary terms vanish (k>=1 when e=1)
G_ = sp.Symbol('g')
Phi4_g = sp.expand(Phi4_g.subs(gg1, G_))
gam = {}                                    # gamma_ab = 1/(a! b!) int s^a g^b
def gm(a, b): return sp.Symbol('gam_%d%d' % (a, b))
total = 0
for term in sp.Add.make_args(Phi4_g):
    c, rest = term.as_independent(s, G_, gp)
    pw = rest.as_powers_dict()
    k, j, e = pw.get(s, 0), pw.get(G_, 0), pw.get(gp, 0)
    if e == 1:                              # int s^k g^j g' = -(k/(j+1)) int s^{k-1} g^{j+1}
        assert k >= 1
        c, k, j = -c*sp.Rational(k, j+1), k-1, j+1
    else:
        assert e == 0
    total += c*sp.factorial(k)*sp.factorial(j)*gm(k, j)
total = sp.expand(total)
print("int_0^inf Phi4 ds =", total)
print("c5 = 8 * that")
# --- step 4: Lemma 5 relations at total degree 4
lam, mu, ep = sp.symbols('lam mu ep')
z = {k: sp.Symbol('z%d' % k) for k in range(2, 6)}
def lgG1(x, order=6):                       # log Gamma(1+x) minus (-euler x): euler cancels in our combination
    return sum((-1)**k*z[k]*x**k/k for k in range(2, order))
# Gamma(-l/2)Gamma(l/2-mu)/Gamma(-mu) = [(-mu)/((-l/2)(l/2-mu))] * exp(lg(-l/2)+lg(l/2-mu)-lg(-mu))
pref = (-mu)/((-lam/2)*(lam/2-mu))
ex = lgG1(-lam/2*ep) + lgG1((lam/2-mu)*ep) - lgG1(-mu*ep)
ser = sp.series(sp.exp(ex), ep, 0, 7).removeO()
# RHS(ep*lam, ep*mu) = 1/2 [ pref*ser/ep + 2/(ep*lam) + 2/(ep*(2mu-lam)) ]; degree-d part = coefficient of ep^d
body = sp.expand(ser)
rhs4 = sp.Rational(1, 2)*pref*body.coeff(ep, 5)       # poles only affect degree -1
rhs4 = sp.factor(sp.simplify(rhs4))
print("RHS degree 4:", sp.expand(rhs4))
# sanity: degree -1 part must vanish: pref*1 + 2/lam + 2/(2mu-lam) = 0
print("pole check:", sp.simplify(pref + 2/lam + 2/(2*mu-lam)))
unk = [gm(a, 4-a) for a in range(0, 4)]    # gamma_{a,4-a}, b>=1
Gdeg = sum(gm(a, 4-a)*(lam**a*mu**(4-a) + (2*mu-lam)**a*mu**(4-a)) for a in range(0, 4))
eqs = sp.Poly(sp.expand(Gdeg - sp.expand(sp.simplify(rhs4))), lam, mu).coeffs()
sol = sp.solve(eqs, unk, dict=True)
print("solution of reflection relations:", sol)
open("sm5_derive_total.txt", "w").write(str(total) + "\n")
if sol:
    sub = sol[0]
    val = sp.expand(total.subs(sub))
    print("int Phi4 (reduced) =", val)
    print("c5 =", sp.expand(8*val))

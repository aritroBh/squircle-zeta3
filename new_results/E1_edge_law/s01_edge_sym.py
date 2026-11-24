# s01: exact symbolic d^3 coefficient of the edge integrand in the RADIAL-PROJECTION parametrisation.
#
# Edge e of K_1 from vertex V0 to V1, y(t) = (1-t)V0 + t V1, t in [0,1].
# x_p(t) = y(t) F_p(y(t))^{-1/p}, F_p(y) = sum_i w_i e^{d lam_i} |<n_i,y>|^p, d = p-1.
# a_i(t) = w_i |<n_i,y(t)>| (affine, sum_i a_i = 1 on the edge), l_i = log|<n_i,y>| + lam_i.
# Exact: |x_p'(t)| = F^{-1/p} sqrt((L - sig*Phi)^2 + h^2 Phi^2), Phi = sum a_i' e^{d l_i} / F,
#   L = edge length, h = distance of the edge line from 0, sig(t) = s0 + L t (coordinate along the edge).
# Endpoint corners: a_0 = al*t (vanishes at V0), a_1 = be*(1-t) (vanishes at V1); all other normals lumped:
#   tau_k = sum_{i in N} a_i(t) l_i^k (tau_0 = 1 - al t - be (1-t)),  gam_k = sum_{i in N} a_i' l_i^k.
import sympy as sp, pickle

t, al, be, L, h, s0, d = sp.symbols('t alpha beta L h s0 d', real=True)
l0, l1 = sp.symbols('l0 l1')
tau = sp.symbols('tau1:4'); gam = sp.symbols('gam1:4')
sig = s0 + L*t
N = 4  # keep d^0..d^3

def mul(a, b):
    return [sp.expand(sum(a[i]*b[k-i] for i in range(k+1))) for k in range(N)]
def inv(a):  # a[0] = 1
    r = [sp.Integer(1)] + [0]*(N-1)
    for k in range(1, N):
        r[k] = sp.expand(-sum(a[i]*r[k-i] for i in range(1, k+1)))
    return r
def log1(a):  # log of series with a[0]=1
    x = [0] + a[1:]
    res = [0]*N; pw = [sp.Integer(1)] + [0]*(N-1)
    for j in range(1, N):
        pw = mul(pw, x)
        res = [sp.expand(res[k] + sp.Rational((-1)**(j+1), j)*pw[k]) for k in range(N)]
    return res
def exp0(a):  # exp of series with a[0]=0
    res = [sp.Integer(1)] + [0]*(N-1); pw = [sp.Integer(1)] + [0]*(N-1)
    for j in range(1, N):
        pw = mul(pw, a)
        res = [sp.expand(res[k] + pw[k]/sp.factorial(j)) for k in range(N)]
    return res

T = [sp.Integer(1)] + [al*t*l0**k + be*(1-t)*l1**k + tau[k-1] for k in range(1, N)]
M = [sp.Integer(0)] + [al*l0**k - be*l1**k + gam[k-1] for k in range(1, N)]
F = [sp.expand(T[k]/sp.factorial(k)) for k in range(N)]
Num = [sp.expand(M[k]/sp.factorial(k)) for k in range(N)]
Phi = mul(Num, inv(F))
oneoverp = [sp.Integer((-1)**k) for k in range(N)]          # 1/(1+d)
Fpow = exp0([sp.expand(-c) for c in mul(oneoverp, log1(F))])  # F^{-1/p}
# sqrt((L - sig Phi)^2 + h^2 Phi^2) = L*sqrt(1+X), X = -2 sig Phi/L + (sig^2+h^2) Phi^2/L^2, Phi = O(d)
Phi2 = mul(Phi, Phi)
X = [sp.expand(-2*sig*Phi[k]/L + (sig**2 + h**2)*Phi2[k]/L**2) for k in range(N)]
sq = [sp.Integer(1)] + [0]*(N-1); pw = [sp.Integer(1)] + [0]*(N-1)
for j in range(1, N):
    pw = mul(pw, X)
    sq = [sp.expand(sq[k] + sp.binomial(sp.Rational(1, 2), j)*pw[k]) for k in range(N)]
S = [sp.expand(L*c) for c in sq]
E = mul(Fpow, S)
E = [sp.expand(c) for c in E]
print("E0 =", sp.simplify(E[0]))
print("E1 =", sp.factor(E[1]))
P3 = sp.Poly(E[3], l0, l1, *tau, *gam)
print("total degree of E3 in logs:", P3.total_degree())
g21 = sp.factor(E[3].coeff(l0, 2).coeff(l1, 1).subs({x: 0 for x in tau + gam}))
g12 = sp.factor(E[3].coeff(l0, 1).coeff(l1, 2).subs({x: 0 for x in tau + gam}))
# confirm no tau/gam in the cross coefficients
print("tau/gam in g21:", E[3].coeff(l0, 2).coeff(l1, 1).free_symbols & set(tau + gam))
print("tau/gam in g12:", E[3].coeff(l0, 1).coeff(l1, 2).free_symbols & set(tau + gam))
print("g21 =", g21)
print("g12 =", g12)
print("g21 polynomial in t of degree", sp.Poly(sp.expand(g21*L**2), t).degree())
I21 = sp.factor(sp.integrate(sp.expand(g21), (t, 0, 1)))
I12 = sp.factor(sp.integrate(sp.expand(g12), (t, 0, 1)))
Itot = sp.factor(I21 + I12)
print("int g21 =", I21)
print("int g12 =", I12)
print("int (g21+g12) =", Itot)
pickle.dump(dict(E=E, g21=g21, g12=g12, I21=I21, I12=I12, Itot=Itot), open('s01_E.pkl', 'wb'))
with open('s01_out.txt', 'w') as f:
    f.write(f"g21 = {g21}\ng12 = {g12}\nint g21 = {I21}\nint g12 = {I12}\nint(g21+g12) = {Itot}\n")
print("saved s01_E.pkl, s01_out.txt")

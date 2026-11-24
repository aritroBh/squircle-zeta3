# step_01.py: extend E1 s01_edge_sym.py to N=6 (d^0..d^5).
# REDUCED run (tau=gam=0): exact cross coefficients g_ab(t) = [l0^a l1^b] E5, a+b=5,
# which by the structure argument (each tau/gam factor consumes one of the 5 d-slots, so it
# cannot appear in a degree-5 l0/l1 monomial) carry ALL the zeta(5) and zeta(2)zeta(3) content.
import sympy as sp, pickle, time

t, al, be, L, h, s0 = sp.symbols('t alpha beta L h s0', real=True)
l0, l1 = sp.symbols('l0 l1')
N = 6

def mul(a, b):
    return [sp.expand(sum(a[i]*b[k-i] for i in range(k+1))) for k in range(N)]
def inv(a):
    r = [sp.Integer(1)] + [0]*(N-1)
    for k in range(1, N):
        r[k] = sp.expand(-sum(a[i]*r[k-i] for i in range(1, k+1)))
    return r
def log1(a):
    x = [0] + a[1:]
    res = [0]*N; pw = [sp.Integer(1)] + [0]*(N-1)
    for j in range(1, N):
        pw = mul(pw, x)
        res = [sp.expand(res[k] + sp.Rational((-1)**(j+1), j)*pw[k]) for k in range(N)]
    return res
def exp0(a):
    res = [sp.Integer(1)] + [0]*(N-1); pw = [sp.Integer(1)] + [0]*(N-1)
    for j in range(1, N):
        pw = mul(pw, a)
        res = [sp.expand(res[k] + pw[k]/sp.factorial(j)) for k in range(N)]
    return res

t0 = time.time()
T = [sp.Integer(1)] + [al*t*l0**k + be*(1-t)*l1**k for k in range(1, N)]
M = [sp.Integer(0)] + [al*l0**k - be*l1**k for k in range(1, N)]
F = [sp.expand(T[k]/sp.factorial(k)) for k in range(N)]
Num = [sp.expand(M[k]/sp.factorial(k)) for k in range(N)]
Phi = mul(Num, inv(F))
oneoverp = [sp.Integer((-1)**k) for k in range(N)]
Fpow = exp0([sp.expand(-c) for c in mul(oneoverp, log1(F))])
Phi2 = mul(Phi, Phi)
sig = s0 + L*t
X = [sp.expand(-2*sig*Phi[k]/L + (sig**2 + h**2)*Phi2[k]/L**2) for k in range(N)]
sq = [sp.Integer(1)] + [0]*(N-1); pw = [sp.Integer(1)] + [0]*(N-1)
for j in range(1, N):
    pw = mul(pw, X)
    sq = [sp.expand(sq[k] + sp.binomial(sp.Rational(1, 2), j)*pw[k]) for k in range(N)]
S = [sp.expand(L*c) for c in sq]
E = mul(Fpow, S)
E = [sp.expand(c) for c in E]
print("E0 =", sp.simplify(E[0]), " E1 =", sp.factor(E[1]))
E5 = E[5]
print("time to build E5 (reduced): %.1f s, terms: %d" % (time.time()-t0, len(sp.Add.make_args(E5))))
# total log degree check
P5 = sp.Poly(E5, l0, l1)
print("total degree of reduced E5 in (l0,l1):", P5.total_degree())
g = {}
for a in range(1, 5):
    b = 5 - a
    g[(a, b)] = sp.expand(E5.coeff(l0, a).coeff(l1, b))
    print(f"g{a}{b} is poly in t of degree", sp.Poly(g[(a, b)], t).degree(), ", terms:", len(sp.Add.make_args(g[(a, b)])))
# integrals over t (polynomial integrals, exact)
Ig = {k: sp.expand(sp.integrate(v, (t, 0, 1))) for k, v in g.items()}
pickle.dump(dict(E=[str(c) for c in E], g={str(k): str(v) for k, v in g.items()},
                 Ig={str(k): str(v) for k, v in Ig.items()}), open('steps/step_01_reduced.pkl', 'wb'))
with open('steps/step_01_out.txt', 'w') as f:
    for k in sorted(g):
        f.write(f"g{k} = {g[k]}\nint g{k} = {Ig[k]}\n")
print("saved steps/step_01_reduced.pkl, steps/step_01_out.txt")

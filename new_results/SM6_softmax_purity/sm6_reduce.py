# SM6 item 2: exact reduction of c6, c7 via gamma_ab and Lemma 5 (sympy).
import sympy as sp
lam, mu, k = sp.symbols('lambda mu k')
D = 6
p = sp.Symbol('p')   # pi^2
Z = {n: sp.Symbol('z%d' % n) for n in range(2, D+3)}   # odd zetas independent symbols
Z[2] = p/6; Z[4] = p**2/90; Z[6] = p**3/945
# unknowns gamma_ab, b>=1, a+b<=D
g = {(a, b): sp.Symbol('g%d%d' % (a, b)) for a in range(D) for b in range(1, D+1-a)}
G = lambda l, m: sum(g[a, b]*l**a*m**b for (a, b) in g)
# RHS = 2mu/(lam(lam-2mu)) (Q-1), logQ = sum_k zeta(k)/k [(lam/2)^k+(mu-lam/2)^k-mu^k]
eps = sp.Symbol('eps')
logQ = sum(Z[n]/n*((lam/2)**n+(mu-lam/2)**n-mu**n) for n in range(2, D+2))
# expand exp(logQ)-1 to total degree D+1 (since dividing by degree-2 lam*(lam-2mu) and times mu -> degree D)
Qm1 = sp.expand(sp.series(sp.exp(eps*0+1*sp.Symbol('T')), sp.Symbol('T'), 0, 1).removeO()*0)  # dummy
T = sp.Symbol('T')
L = sp.expand(logQ.subs({lam: T*lam, mu: T*mu}))
E = sp.series(sp.exp(L), T, 0, D+2).removeO() - 1
E = sp.expand(E)
Rhs = 0
for deg in range(2, D+2):
    c = sp.expand(E.coeff(T, deg))
    q, r = sp.div(sp.Poly(c, lam, mu), sp.Poly(lam*(lam-2*mu), lam, mu))
    assert r.is_zero, (deg, r)
    Rhs += sp.expand(2*mu*q.as_expr())      # degree deg-2+1 = deg-1
# coefficient comparison: G(lam,mu)+G(2mu-lam,mu) - Rhs = 0 up to total degree D
expr = sp.expand(G(lam, mu)+G(2*mu-lam, mu)-Rhs)
P = sp.Poly(expr, lam, mu)
eqs = [c for (a, b), c in P.terms() if a+b <= D for c in [c]]
eqs = [c for (m_, c) in [((a, b), c) for (a, b), c in P.terms()] if True]
eqs = [c for (ab, c) in P.terms() if sum(ab) <= D]
# elementary gamma_a1 = eta(a+2)/2^(a+1), eta(n)=(1-2^(1-n)) zeta(n)
for a in range(D):
    eqs.append(g[a, 1]-(1-sp.Rational(2)**(1-(a+2)))*Z[a+2]/sp.Rational(2)**(a+1))
for d in range(1, D+1):
    e2 = [c for (ab, c) in P.terms() if sum(ab) == d] + [g[a, 1]-(1-sp.Rational(2)**(-1-a))*Z[a+2]/sp.Rational(2)**(a+1) for a in range(d)]
    r = sp.solve(e2, [g[a, b] for (a, b) in g if a+b <= d], dict=True)
    print(d, len(r))
    if d==3: print(e2, r)
sol = sp.solve(eqs, list(g.values()), dict=True)
assert len(sol) == 1
sol = sol[0]
free = [v for v in g.values() if v not in sol]
print("free gamma-combinations (solve returns unsolved symbols):", free)
import pickle
pickle.dump((sol, free), open("sm6_reduce.pkl", "wb"))
def red(e): return sp.expand(e.subs(sol))
c = {}
c[4] = 8*(g[0, 3]+g[1, 2])
c[5] = 8*(g[0, 4]+g[1, 3]+3*g[3, 1])
c[6] = 8*(g[0, 5]+g[1, 4]+3*g[3, 2])
c[7] = 8*(g[0, 6]+g[1, 5]+3*g[3, 3]+45*g[5, 1])
for n in (4, 5, 6, 7):
    print("c%d = %s" % (n, c[n]))
    print("   reduced: ", sp.factor(red(c[n])) if False else red(c[n]))
# which free symbols appear
for n in (4, 5, 6, 7):
    print("c%d free symbols after reduction:" % n, red(c[n]).free_symbols & set(free))
# the general solution in terms of free
print("solved:")
for a in range(D):
    for b in range(1, D+1-a):
        if g[a, b] in sol: print(g[a, b], "=", sol[g[a, b]])

"""
compact_forms.py -- extracts theta_{ij} (Taylor coefficients of Theta(W,mu)) and the tail constants T_k,
and verifies the compact representation
      Q_k(Xi) = sum_{i=1}^k theta_{i,k-i} M_i(Xi) + T_k,   M_i(Xi) = E[(Xi + Y)^i],  Y = xi - k1,
against the brute-force symbolic expansion of symbolic_Qk.py (exact sympy identity, k <= 4).
"""
import sympy as sp
import symbolic_Qk as S

mu, Xi, W, h = S.mu, S.Xi, S.W, S.h
a, b, c, d = S.a, S.b, S.c, S.d
k2, k3, k4 = S.k2, S.k3, S.k4
K = S.K
Tfull = sp.expand(S.T.subs(h, 1))
P = sp.Poly(Tfull, W, mu)
theta = {}
for (i, j), t in P.terms():
    theta[(i, j)] = sp.factor(t)
cm = {0: 1, 1: 0, 2: k2, 3: k3, 4: k4 + 3*k2**2}
def M(i):
    return sp.expand(sum(sp.binomial(i, r)*Xi**(i - r)*cm[r] for r in range(i + 1)))
tail = sp.expand(-S.tail)
print('theta_{ij} (coefficient of W^i mu^j in Theta):')
for k in range(1, K + 1):
    for i in range(k, -1, -1):
        j = k - i
        if (i, j) in theta:
            print('   theta_%d%d =' % (i, j), theta[(i, j)])
print('tail constants T_k = [mu^k](-Tail):')
Tk = {}
for k in range(1, K + 1):
    Tk[k] = sp.factor(tail.coeff(mu, k))
    print('   T_%d =' % k, Tk[k])
print('check Q_k == sum_i theta_{i,k-i} M_i(Xi) + T_k :')
for k in range(1, K + 1):
    comp = sum(theta.get((i, k - i), 0)*M(i) for i in range(1, k + 1)) + Tk[k]
    print('   k=%d  difference =' % k, sp.simplify(sp.expand(comp - S.Q[k])))

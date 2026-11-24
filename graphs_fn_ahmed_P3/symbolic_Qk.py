"""
symbolic_Qk.py -- exact symbolic expansion coefficients Q_k(Xi) for general f.

  D_n ~ sum_{k>=1} mu^k Q_k(Xi),   mu = 1/(n-1),   Xi = log(e f'(1) n / 2) = log n + log a + 1 - log 2.

Q_k is a polynomial of degree k in Xi whose coefficients are rational functions of
a=f'(1), b=f''(1), c=f'''(1), d=f''''(1), ... and the cumulants
  k2 = 1 + pi^2/6,  k3 = 2 - 3 zeta(3)/2,  k4 = 6 + pi^4/15,  ...
of xi = log(1/S), S ~ g'(s) ds.

Derivation (see result.md, Theorem 3):
  main part  = E[ Ttil(mu Xi + mu (xi - kappa), mu) ]            (Ttil solves A(T) - mu Btil(T) = W)
  tail part  = sum_{i,j} ttil_ij mu^{i+j} (-1)^i i! sum_k b_k (n a)^{1-2k} / (2k-1)^{i+1},
               b_k = -binom(-1/2,k)  (g'(s) = sum_k b_k s^{-2k}),  1/(n a) = mu/(a (1+mu)).
Outputs: LaTeX for Q_1..Q_4, checks f=x against Ahmed (exact formula), checks f=x^m.
"""
import sympy as sp

K = 4
a, b, c, d, e5 = sp.symbols('a b c d e5', positive=True)
jet = [None, a, b, c, d, e5]
tau, W, mu, Xi = sp.symbols('tau W mu Xi')
k2, k3, k4, k5 = sp.symbols('k2 k3 k4 k5')

# f(1-tau) = sum_j (-1)^j f^{(j)}(1) tau^j / j!
f1t = 1 + sum((-1)**j*jet[j]*tau**j/sp.factorial(j) for j in range(1, K + 2))
fp1t = sum((-1)**(j - 1)*jet[j]*tau**(j - 1)/sp.factorial(j - 1) for j in range(1, K + 2))
A = sp.series(-sp.log(f1t), tau, 0, K + 1).removeO()
Btil = sp.series(sp.log(fp1t/a), tau, 0, K + 1).removeO()

# solve A(T) - mu*Btil(T) = W by iteration with a bookkeeping parameter h (W->hW, mu->h mu)
h = sp.Symbol('h')
T = 0
for it in range(K + 1):
    AT = sp.expand(A.subs(tau, T))
    BT = sp.expand(Btil.subs(tau, T))
    a1 = a
    T = sp.expand((h*W + h*mu*BT - (AT - a1*T))/a1)
    T = sp.series(T, h, 0, K + 1).removeO()
T = sp.expand(T)

# main part: E[T(mu Xi + mu Y, mu)], Y = xi - kappa centered, E[Y]=0, E[Y^2]=k2, ...
# central moments from cumulants
cm = {0: 1, 1: 0, 2: k2, 3: k3, 4: k4 + 3*k2**2, 5: k5 + 10*k3*k2}
Y = sp.Symbol('Y')
Tmain = sp.expand(T.subs(h, 1).subs(W, mu*Xi + mu*Y))
Tmain = sum(Tmain.coeff(mu, k)*mu**k for k in range(0, K + 1))
poly = sp.Poly(Tmain, Y)
main = 0
for (p,), coeff in poly.terms():
    main += coeff*cm[p]
main = sp.expand(main)

# tail part
Tpoly = sp.Poly(sp.expand(T.subs(h, 1)), W, mu)
tail = 0
for (i, j), t in Tpoly.terms():
    k0 = i + j
    for kk in range(1, K + 1):
        p = 2*kk - 1
        if k0 + p > K:
            break
        bk = -sp.binomial(sp.Rational(-1, 2), kk)
        Na_inv = mu/(a*(1 + mu))
        tail += t*mu**k0*(-1)**i*sp.factorial(i)*bk*Na_inv**p/sp.Integer(2*kk - 1)**(i + 1)
tail = sp.series(tail, mu, 0, K + 1).removeO()

D = sp.expand(main - tail)
Q = {}
for k in range(1, K + 1):
    Q[k] = sp.collect(sp.expand(D.coeff(mu, k)), Xi)

def show(Qk):
    out = []
    pk = sp.Poly(Qk, Xi)
    for (m,), coeff in sorted(pk.terms(), reverse=True):
        out.append((m, sp.factor(sp.simplify(coeff))))
    return out

if __name__ == '__main__':
    for k in range(1, K + 1):
        print('Q_%d:' % k)
        for m, coeff in show(Q[k]):
            print('   Xi^%d :' % m, coeff)
    # ---- check f = x (a=1, b=c=d=0) against the mu-expansion of Ahmed's exact formula
    nn = 1 + 1/mu
    kap = 1 - sp.log(2)
    # log R(eps) = sum_j a_j eps^j ; k_j = j! a_j (-1)^j
    L = sp.Symbol('L')
    epsv = sp.Symbol('eps')
    lnR = (sp.log(2) - 1)*epsv + sum((((sp.zeta(j) + 1) if j % 2 == 0 else (sp.zeta(j)*(1 - sp.Integer(2)**(1 - j)) - 1))/j)*epsv**j for j in range(2, K + 2))
    kj = {j: sp.factorial(j)*lnR.coeff(epsv, j)*(-1)**j for j in range(1, K + 2)}
    print('cumulants:', [sp.nsimplify(kj[j]) for j in range(1, K + 1)])
    vals = {k2: kj[2], k3: kj[3], k4: kj[4], k5: kj[5]}
    # Ahmed: D = 1 - n^{-mu} R(mu) + sum_k binom(1/2,k)/(n^{2k-1}((2k-1)(n-1)-1))
    Dah = 1 - sp.exp(-mu*L)*sp.exp(lnR.subs(epsv, mu))
    Dah += sum(sp.binomial(sp.Rational(1, 2), kk)/(nn**(2*kk - 1)*((2*kk - 1)/mu - 1)) for kk in range(1, 3))
    # expand in mu with L treated as O(1)
    Dah_s = sp.series(Dah, mu, 0, K + 1).removeO()
    ours = sp.expand(D.subs({a: 1, b: 0, c: 0, d: 0, e5: 0}).subs(vals).subs(Xi, L + kap))
    for k in range(1, K + 1):
        diff = sp.simplify(sp.expand(Dah_s.coeff(mu, k) - ours.coeff(mu, k)))
        print('f=x, order mu^%d: ours - Ahmed =' % k, diff)
    # ---- check f = x^m  (D_f(n) = D_x(m n)): use a=m, b=m(m-1), c=m(m-1)(m-2), ...
    msym = sp.Symbol('m', positive=True)
    sub_m = {a: msym, b: msym*(msym - 1), c: msym*(msym - 1)*(msym - 2), d: msym*(msym - 1)*(msym - 2)*(msym - 3), e5: 0}
    mu2 = sp.Symbol('mu2')  # 1/(m n - 1) for x
    # Ahmed for N = m n, eps = 1/(N-1);  N - 1 = m(n-1) + m - 1  -> eps = mu/(m + (m-1) mu)
    epsN = mu/(msym + (msym - 1)*mu)
    NN = msym*nn
    DahN = 1 - sp.exp(-epsN*(L + sp.log(msym)))*sp.exp(lnR.subs(epsv, epsN))
    DahN += sum(sp.binomial(sp.Rational(1, 2), kk)/(NN**(2*kk - 1)*((2*kk - 1)/epsN - 1)) for kk in range(1, 3))
    DahN_s = sp.series(DahN, mu, 0, K + 1).removeO()
    oursm = sp.expand(D.subs(sub_m).subs(vals).subs(Xi, L + sp.log(msym) + kap))
    for k in range(1, K + 1):
        diff = sp.simplify(sp.expand(DahN_s.coeff(mu, k) - oursm.coeff(mu, k)))
        print('f=x^m, order mu^%d: ours - Ahmed(mn) =' % k, diff)
    # ---- LaTeX
    print()
    for k in range(1, K + 1):
        print('Q_%d LaTeX:' % k, sp.latex(sp.collect(sp.factor_terms(Q[k]), Xi)))

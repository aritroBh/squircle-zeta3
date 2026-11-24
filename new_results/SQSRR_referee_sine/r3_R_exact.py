# Part 2: exact evaluation of R = sum_{n>=3} b_n'(2) A_n, and of sum b_n'C_n, sum b_n'^2, from the law C1.
# Pipeline (every step an identity; see REFEREE.md Sec. 4):
#  (A) n = 4k+eps, k>=1:  beta_n A_n = r1(k) + rH(k) H_k + rG(k) H_{2k} + rO(k) O_k,  O_k = sum_{i<=k} 1/(2i-1)^2
#  (B) partial fractions in k: poles only at k = 0, -eps/2
#  (C) each term c X_k / D^j, D in {k, 2k+1, 2k-1}, rewritten as a parity-weighted sum over m (m = 2k, 2k+1, 2k-1)
#      of Y_m / m^j, Y in {1, H_m, A_m, H2_m, A2_m}  (A_m = sum_{l<=m} (-1)^{l+1}/l, A2 likewise), plus rational
#      sums with poles at m=0,-1.  All truncations are at M = 2N + O(1); T := H_M (formal, must cancel).
#  (D) P(Y,s,j) = sum_m s^m Y_m/m^j  ->  depth-2 sums Z and depth-1 sums L
#  (E) Z reduced by dshuffle.reduce_weight (double shuffle + distribution; one free weight-4 element Z(3-,1-)).
import sympy as sp, pickle
from dshuffle import reduce_weight, L, ZS, l2, zsym, T
k, m = sp.symbols('k m', positive=True)
EVEN, ODD = 'even', 'odd'
def Lreg(a, s):
    return T if (a == 1 and s == 1) else L(a, s)
def Zreg(a, s, b, t):
    if not (a == 1 and s == 1): return zsym(a, s, b, t)
    if b == 1 and t == 1: return (T**2 - ZS[2])/2
    return T*L(b, t) - zsym(b, t, 1, 1) - L(b+1, t)
def P(Y, s, j):
    if Y == '1':  return Lreg(j, s)
    if Y == 'H':  return Zreg(j, s, 1, 1) + Lreg(j+1, s)
    if Y == 'A':  return -Zreg(j, s, 1, -1) - Lreg(j+1, -s)
    if Y == 'H2': return Zreg(j, s, 2, 1) + Lreg(j+2, s)
    if Y == 'A2': return -Zreg(j, s, 2, -1) - Lreg(j+2, -s)
def Spar(par, Y, j):      # sum_{m>=1} psi_par(m) Y_m / m^j
    sg = 1 if par == EVEN else -1
    return sp.Rational(1, 2)*(P(Y, 1, j) + sg*P(Y, -1, j))
def Srat(par, f):         # sum_{m>=1} psi_par(m) f(m), f rational with poles at 0, -1
    out = 0
    for term in sp.Add.make_args(sp.apart(sp.together(f), m)):
        num, den = sp.fraction(sp.factor(term))
        pw = sp.Poly(den, m)
        roots = sp.roots(pw)
        assert len(roots) == 1, term
        (r0, j), = roots.items()
        c = sp.simplify(term*(m-r0)**j); assert not c.has(m)
        if r0 == 0: out += c*Spar(par, '1', j)
        else:
            assert r0 == -1
            flip = ODD if par == EVEN else EVEN     # m' = m+1 >= 2
            out += c*(Spar(flip, '1', j) - (1 if flip == ODD else 0))
    return out
# map X_k for each D -> list of (Y, coeff) + rational extra g(m)
def Xmap(X, D):
    if D == 'k':
        return {'1': [('1', 1)], 'Hk': [('H', 1), ('A', -1)], 'H2k': [('H', 1)],
                'O': [('H2', sp.Rational(1, 2)), ('A2', sp.Rational(1, 2))]}[X], 0
    if D == '2k+1':
        return {'1': [('1', 1)], 'Hk': [('H', 1), ('A', -1)], 'H2k': [('H', 1)],
                'O': [('H2', sp.Rational(1, 2)), ('A2', sp.Rational(1, 2))]}[X], \
               {'1': 0, 'Hk': 0, 'H2k': -1/m, 'O': -1/m**2}[X]
    if D == '2k-1':
        return {'1': [('1', 1)], 'Hk': [('H', 1), ('A', -1)], 'H2k': [('H', 1)],
                'O': [('H2', sp.Rational(1, 2)), ('A2', sp.Rational(1, 2))]}[X], \
               {'1': 0, 'Hk': 2/(m+1), 'H2k': 1/(m+1), 'O': 0}[X]
def sum_k(r, X):
    """exact (regularised) value of sum_{k>=1} r(k) X_k"""
    out = 0
    for term in sp.Add.make_args(sp.apart(sp.together(r), k)):
        num, den = sp.fraction(sp.factor(term))
        roots = sp.roots(sp.Poly(den, k)); assert len(roots) == 1, term
        (r0, j), = roots.items()
        c = sp.simplify(term*(k-r0)**j); assert not c.has(k)
        if r0 == 0:   D, cm = 'k', c*2**j                                   # 1/k^j = 2^j/m^j, m=2k
        elif r0 == sp.Rational(-1, 2): D, cm = '2k+1', c*2**j               # 1/(k+1/2)^j = 2^j/m^j
        elif r0 == sp.Rational(1, 2):  D, cm = '2k-1', c*2**j
        else: raise ValueError(term)
        par = EVEN if D == 'k' else ODD
        ys, g = Xmap(X, D)
        for Y, cy in ys: out += cm*cy*Spar(par, Y, j)
        if g != 0: out += Srat(par, cm*g/m**j)
        if D == '2k+1' and X == '1': out -= cm     # remove m = 1 (k = 0) term
    return sp.expand(out)
def classes():
    for eps in (1, -1):
        n = 4*k+eps; a = 4*k; c = 4*k+2*eps
        beta = sp.Rational(2)*eps/(a*c**2)
        rho = -2*n*(3*n-1)/(a**2*c**3) if eps == 1 else 2*(3*n**2+3*n-2)/(a**2*c**3)
        C = -sp.Rational(eps)*n/(4*a**2*c)
        yield eps, n, a, c, beta, rho, C
def evaluate(expr):
    """substitute depth-2 reductions"""
    expr = sp.expand(expr)
    for w in range(2, 6):
        red, free = reduce_weight(w)
        expr = sp.expand(expr.subs(red))
    return sp.expand(expr.subs(ZS[2]**2, sp.Rational(5, 2)*ZS[4]))
if __name__ == '__main__':
    R = 0; SC = 0; SQ = 0
    for eps, n, a, c, beta, rho, C in classes():
        R += sum_k(beta*rho, '1') + sum_k(beta*(-2*eps/a**3), 'Hk') \
             + sum_k(beta*2*eps*(3*n**2+2*eps*n+3)/(a**3*c**2), 'H2k') + sum_k(beta*2*eps*n/(a**2*c), 'O')
        SC += sum_k(beta*C, '1'); SQ += sum_k(beta**2, '1')
    out = {}
    for name, v in (('sum beta C', SC), ('sum_{n>=3} beta^2', SQ), ('R', R)):
        v = evaluate(v)
        print(name, ': T-dependence:', sp.simplify(v.coeff(T, 1)), sp.simplify(v.coeff(T, 2)))
        v = sp.expand(v.subs(T, 0))
        print('   =', v)
        out[name] = v
    pickle.dump(out, open('r3_R_exact.pkl', 'wb'))
    # compare with C3
    pi2 = 6*ZS[2]; pi4 = 90*ZS[4]
    C3 = sp.Rational(5, 16) - sp.Rational(13, 96)*pi2 + sp.Rational(47, 64)*ZS[3] + pi2*l2/32 - pi4/384 \
         + sp.Rational(31, 256)*ZS[5] + pi2*ZS[3]/256
    print('R - C3 =', sp.simplify(sp.expand(out['R'] - C3)))
    print('sum beta C - (1-zeta3)/64 =', sp.simplify(out['sum beta C'] - (1-ZS[3])/64))

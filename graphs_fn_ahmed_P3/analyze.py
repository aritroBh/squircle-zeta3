"""
analyze.py -- compares the high-precision defects D_n (data_defects.json) with the asymptotic
expansion of Theorem 3, extracts the constant C(f) = log f'(1) + 1 - log 2, and performs the
Richardson-type (linear least-squares-free, exactly determined) extrapolations of the remainder.

Prints:
  [0] quadrature validation vs Ahmed's exact Gamma-function formula (f=x, and f=x^2 via D_{x^2}(n)=D_x(2n))
  [1] closed-form Q_2, Q_3 (symbolic_Qk.py) vs the numeric engine (fnlib.expansion_polys)
  [2] ratio test  r_K(n)/(mu^{K+1} P_{K+1}(log n)) -> 1   (r_K = D_n - sum_{k<=K} mu^k P_k)
  [3] constant extraction from single n using closed-form Q_2,Q_3 (and engine to order 12)
  [4] Richardson A: Q_2,Q_3 fixed (proved), C and all P_4..P_7 coefficients unknown, solved exactly
  [5] Richardson B: theory-free fit, P_1..P_7 all unknown (35 unknowns)
  [6] slope-2/e principle: L_n = 1 + x_n^* + c_2(f) mu^2 + O(mu^3 log n)
"""
import json
import mpmath as mp
import sympy as sp
from fnlib import TESTS, make_callables, jets, expansion_polys, eval_expansion, R, moments_xi

mp.mp.dps = 130
data = json.load(open('data_defects.json'))
NAMES = ['x', 'x^2', 'expm1', 'x+x^3', 'mobius', 'sin4']
xs = sp.Symbol('x', positive=True)
KAP = 1 - mp.log(2)
K2 = 1 + mp.pi**2/6
K3 = 2 - mp.mpf(3)/2*mp.zeta(3)

def Dtab(name):
    return {int(r['n']): mp.mpf(r['D']) for r in data if r['name'] == name}

def L_ahmed(n):
    n = mp.mpf(n); eps = 1/(n - 1)
    s = mp.nsum(lambda k: mp.binomial(mp.mpf(1)/2, k)/(n**(2*k - 1)*((2*k - 1)*(n - 1) - 1)), [1, mp.inf])
    return 1 + R(eps)/n**eps - s

def Q2(a, b, Xi):
    return -((a**2 - b)/(2*a**3)*Xi**2 + b/a**3*Xi + ((a**2 - b)*K2 - a)/(2*a**3))

def Q3(a, b, c, Xi):
    return ((a**4 - 3*a**2*b - a*c + 3*b**2)/(6*a**5)*Xi**3
            + (3*a**2*b + a*c - 4*b**2)/(2*a**5)*Xi**2
            + (a**4*K2 - 3*a**2*b*K2 - a*c*K2 + 3*b**2*K2 + 2*b**2)/(2*a**5)*Xi
            + (a**4*K3 + 9*a**2*b*K2 - 3*a**2*b*K3 - 6*a*b + 3*a*c*K2 - a*c*K3 - 12*b**2*K2 + 3*b**2*K3)/(6*a**5))

def digits(err, ref):
    if err == 0:
        return 999
    return float(-mp.log10(abs(err/ref)))

KENG = 12
info = {}
for name in NAMES:
    expr = TESTS[name]
    a, b, c = [mp.mpf(sp.N(sp.diff(expr, xs, k).subs(xs, 1), 150)) for k in (1, 2, 3)]
    A, B = jets(expr, KENG + 2)
    P = expansion_polys(A, B, KENG)
    info[name] = dict(a=a, b=b, c=c, P=P, D=Dtab(name), C=mp.log(a) + KAP)

print('=' * 100)
print('[0] quadrature vs Ahmed exact formula (relative residuals)')
worst0 = 0
for n in sorted(info['x']['D']):
    r1 = (info['x']['D'][n] - (2 - L_ahmed(n)))/info['x']['D'][n]
    r2 = (info['x^2']['D'][n] - (2 - L_ahmed(2*n)))/info['x^2']['D'][n]
    worst0 = max(worst0, abs(r1), abs(r2))
print('   max relative residual over n=10^3..10^20, f=x and f=x^2:', mp.nstr(worst0, 3))

print('=' * 100)
print('[1] closed-form Q_2,Q_3 vs numeric engine at Xi = 3.7 (residual)')
for name in NAMES:
    I = info[name]
    a, b, c, P = I['a'], I['b'], I['c'], I['P']
    L = mp.mpf('3.7') - I['C']          # log n with Xi = 3.7
    Xi = mp.mpf('3.7')
    e2 = mp.polyval(P[2][::-1], L) - Q2(a, b, Xi)
    e3 = mp.polyval(P[3][::-1], L) - Q3(a, b, c, Xi)
    print('   %-7s a=%-10s b=%-10s c=%-10s  |dQ2|=%s |dQ3|=%s' % (name, mp.nstr(a, 6), mp.nstr(b, 6), mp.nstr(c, 6),
          mp.nstr(abs(e2), 3), mp.nstr(abs(e3), 3)))

print('=' * 100)
print('[2] ratio test: r_K(n) / (mu^{K+1} P_{K+1}(log n))  (should -> 1, deviation O(mu log n))')
for name in NAMES:
    I = info[name]
    print('  f = %s' % name)
    for n in [10**6, 10**10, 10**15, 10**20]:
        D = I['D'][n]
        mu = 1/mp.mpf(n - 1); L = mp.log(n)
        row = []
        for K in range(1, 9):
            rK = D - eval_expansion(I['P'], n, K)
            nxt = mu**(K + 1)*mp.polyval(I['P'][K + 1][::-1], L)
            row.append(mp.nstr(rK/nxt - 1, 2))
        print('    n=1e%-3d ' % len(str(n)[1:]) + ' '.join('%10s' % s for s in row))

print('=' * 100)
print('[3] single-n constant extraction  C_est(n) = (a/mu)[D_n - sum_{k=2}^K mu^k P_k] - log n ;  exact C = log a + 1 - log 2')
for name in NAMES:
    I = info[name]
    a, b, c = I['a'], I['b'], I['c']
    print('  f = %-7s  C(f) = %s' % (name, mp.nstr(I['C'], 40)))
    for n in [10**6, 10**8, 10**12, 10**16, 10**20]:
        D = I['D'][n]; mu = 1/mp.mpf(n - 1); L = mp.log(n)
        Xi = L + I['C']
        c3 = (a/mu)*(D - mu**2*Q2(a, b, Xi) - mu**3*Q3(a, b, c, Xi)) - L
        c12 = (a/mu)*(D - (eval_expansion(I['P'], n, KENG) - mu*mp.polyval(I['P'][1][::-1], L))) - L
        print('     n=1e%-3d  closed-form Q2,Q3: digits=%6.1f  resid=%-10s | engine K=12: digits=%6.1f resid=%s' % (
            len(str(n)) - 1, digits(c3 - I['C'], I['C']), mp.nstr(c3 - I['C'], 3),
            digits(c12 - I['C'], I['C']), mp.nstr(c12 - I['C'], 3)))
    I['C_single'] = (a/(1/mp.mpf(10**20 - 1)))*(I['D'][10**20] - (eval_expansion(I['P'], 10**20, KENG)
                     - mp.polyval(I['P'][1][::-1], mp.log(10**20))/(10**20 - 1))) - mp.log(10**20)

def fit(name, nset, kmin, kmax, fixed):
    """Solve exactly: D_n - fixed(n) = c0*mu/a + sum_{k=kmin}^{kmax} sum_{i=0}^k c_{k,i} mu^k L^i ;
    returns c0 (the constant C) and the coefficient dict."""
    I = info[name]; a = I['a']
    cols = [('C', None)] + [(k, i) for k in range(kmin, kmax + 1) for i in range(k + 1)]
    assert len(nset) >= len(cols)
    nset = nset[:len(cols)]
    M = mp.matrix(len(cols), len(cols)); rhs = mp.matrix(len(cols), 1)
    for r, n in enumerate(nset):
        mu = 1/mp.mpf(n - 1); L = mp.log(n)
        for cidx, col in enumerate(cols):
            if col[0] == 'C':
                M[r, cidx] = mu/a
            else:
                k, i = col
                M[r, cidx] = mu**k*L**i
        rhs[r] = I['D'][n] - fixed(n, mu, L)
    sol = mp.lu_solve(M, rhs)
    return sol[0], {cols[j]: sol[j] for j in range(1, len(cols))}

print('=' * 100)
print('[4] Richardson A: Q_1 log-coefficient, Q_2, Q_3 fixed (proved); C and P_4..P_7 (22 coeffs) unknown')
allns = sorted(info['x']['D'])
for name in NAMES:
    I = info[name]; a, b, c = I['a'], I['b'], I['c']
    def fixedA(n, mu, L, I=I, a=a, b=b, c=c):
        Xi = L + I['C']   # Q2,Q3 contain C only through Xi (proved closed forms)
        return mu*L/a + mu**2*Q2(a, b, Xi) + mu**3*Q3(a, b, c, Xi)
    ns = [n for n in allns if n >= 10**9]
    Cfit, coef = fit(name, ns, 4, 7, fixedA)
    p4 = I['P'][4]
    dev4 = max(abs(coef[(4, i)] - p4[i])/max(abs(p4[i]), mp.mpf(10)**-30) for i in range(5))
    I['C_richA'] = Cfit
    print('   %-7s C_fit - C = %-10s (digits %.1f);  fitted P_4 coeffs vs engine: max rel dev %s' % (
        name, mp.nstr(Cfit - I['C'], 3), digits(Cfit - I['C'], I['C']), mp.nstr(dev4, 3)))

print('=' * 100)
print('[5] Richardson B (theory-free): D_n = sum_{k=1}^{7} mu^k sum_{i<=k} c_{k,i} L^i, all 35 coefficients unknown')
for name in NAMES:
    I = info[name]; a = I['a']
    ns = [n for n in allns if n >= 10**4]
    cols = [(k, i) for k in range(1, 8) for i in range(k + 1)]
    # use 35 values, spread
    step = max(1, len(ns)//len(cols))
    nsel = ns[-len(cols):]
    M = mp.matrix(len(cols), len(cols)); rhs = mp.matrix(len(cols), 1)
    for r, n in enumerate(nsel):
        mu = 1/mp.mpf(n - 1); L = mp.log(n)
        for j, (k, i) in enumerate(cols):
            M[r, j] = mu**k*L**i
        rhs[r] = I['D'][n]
    sol = mp.lu_solve(M, rhs)
    c10, c11 = sol[0], sol[1]
    Cfree = c10/c11      # C = (coefficient of mu) / (coefficient of mu L)
    I['C_richB'] = Cfree
    print('   %-7s 1/a: fitted %s (digits %.1f);  C = c10/c11: resid %s (digits %.1f)' % (
        name, mp.nstr(c11, 12), digits(c11 - 1/a, 1/a), mp.nstr(Cfree - I['C'], 3), digits(Cfree - I['C'], I['C'])))

print('=' * 100)
print('[6] slope-2/e principle: x_n^* := point where (f^n)\' = 2/e ;  L_n - 1 - x_n^* = c_2 mu^2 + O(mu^3 log n)')
print('    c_2(f) = [ (a^2-b) pi^2/12 + (a^2-a-b)/2 ] / a^3')
for name in NAMES:
    I = info[name]; a, b = I['a'], I['b']
    f, fp = make_callables(TESTS[name])
    c2 = ((a**2 - b)*mp.pi**2/12 + (a**2 - a - b)/2)/a**3
    row = []
    for n in [10**4, 10**6, 10**9, 10**12]:
        nn = mp.mpf(n)
        F = lambda t: (nn - 1)*(-mp.log(f(1 - t))) - mp.log(fp(1 - t)) - mp.log(nn) - KAP
        t0 = (mp.log(nn) + KAP + mp.log(a))/(a*(nn - 1))
        tstar = mp.findroot(F, t0)
        Ln = 2 - I['D'][n]
        mu = 1/(nn - 1)
        resid = Ln - 1 - (1 - tstar) - c2*mu**2
        row.append('n=1e%d: (L-1-x*)/mu^2=%s  resid/(mu^3 log n)=%s' % (len(str(n)) - 1,
                   mp.nstr((Ln - 2 + tstar)/mu**2, 12), mp.nstr(resid/(mu**3*mp.log(nn)), 6)))
    print('   %-7s c_2 = %s' % (name, mp.nstr(c2, 20)))
    for s in row:
        print('        ' + s)

print('=' * 100)
print('SUMMARY TABLE (constant term C(f) = log f\'(1) + 1 - log 2)')
print('%-8s %-45s %-12s %-12s %-12s' % ('f', 'C(f) (exact formula, 40 digits)', 'single n=1e20', 'Richardson A', 'theory-free B'))
for name in NAMES:
    I = info[name]
    print('%-8s %-45s %-12.1f %-12.1f %-12.1f' % (name, mp.nstr(I['C'], 40), digits(I['C_single'] - I['C'], I['C']),
          digits(I['C_richA'] - I['C'], I['C']), digits(I['C_richB'] - I['C'], I['C'])))

print('=' * 100)
print('[7] 1/n-form of the second-order term (Corollary 2):')
print('    D_n = (1+log(an/2))/(an) - (a^2-b)/(2a^3) [log^2(an/2) + pi^2/6]/n^2 + 1/(2a^2 n^2) + O(log^3 n/n^3)')
for name in NAMES:
    I = info[name]; a, b = I['a'], I['b']
    row = []
    for n in [10**6, 10**10, 10**15, 10**20]:
        nn = mp.mpf(n); l = mp.log(a*nn/2)
        approx = (1 + l)/(a*nn) - (a**2 - b)/(2*a**3)*(l**2 + mp.pi**2/6)/nn**2 + 1/(2*a**2*nn**2)
        row.append(mp.nstr((I['D'][n] - approx)*nn**3/mp.log(nn)**3, 6))
    print('   %-7s (D_n - approx) n^3/log^3 n at n=1e6,1e10,1e15,1e20:  %s' % (name, '  '.join(row)))

print('=' * 100)
print('[8] numerical values of Q_2(Xi) = q22 Xi^2 + q21 Xi + q20 (closed form) for the test functions')
for name in NAMES:
    I = info[name]; a, b = I['a'], I['b']
    q22 = -(a**2 - b)/(2*a**3); q21 = -b/a**3; q20 = -((a**2 - b)*K2 - a)/(2*a**3)
    print('   %-7s a=%-22s b=%-22s q22=%-22s q21=%-22s q20=%s' % (name, mp.nstr(a, 18), mp.nstr(b, 18),
          mp.nstr(q22, 18), mp.nstr(q21, 18), mp.nstr(q20, 18)))

print('=' * 100)
print('[9] slope-2/e principle to high precision: L_n - 1 - x_n^*  vs  -(sum_k mu^k P_k(L) - sum_ij t_ij mu^{i+j}(L+k1)^i), K=12')
from fnlib import solve_T
for name in NAMES:
    I = info[name]; a = I['a']
    A, B = jets(TESTS[name], KENG + 2)
    Tser = solve_T(A, B, KENG)
    f, fp = make_callables(TESTS[name])
    row = []
    for n in [10**6, 10**9, 10**12]:
        nn = mp.mpf(n); mu = 1/(nn - 1); L = mp.log(nn)
        F = lambda t: (nn - 1)*(-mp.log(f(1 - t))) - mp.log(fp(1 - t)) - mp.log(nn) - KAP
        tstar = mp.findroot(F, (L + KAP + mp.log(a))/(a*(nn - 1)), tol=mp.mpf(10)**-125)
        lhs = (2 - I['D'][n]) - 1 - (1 - tstar)
        tau_ser = sum(t*(mu*(L + KAP))**i*mu**j for (i, j), t in Tser.items())
        rhs = -(eval_expansion(I['P'], n, KENG) - tau_ser)
        row.append('n=1e%d: rel.resid=%s' % (len(str(n)) - 1, mp.nstr((lhs - rhs)/lhs, 3)))
    print('   %-7s %s' % (name, '   '.join(row)))

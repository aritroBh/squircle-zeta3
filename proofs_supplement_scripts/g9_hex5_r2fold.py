# g9_hex5_r2fold.py -- the edge r = 2 of the hexagon |x|^p+|y|^p+3|x-y|^p <= 4^p (corners (1,0),(2,2),(0,1), weights 1,3,1).
# (1) the 13 harmonic polylogarithms at x = -1 are checked numerically (100 digits, own series + Cohen-Villegas-Zagier acceleration)
#     against their values in the 20 alternating-MZV constants (proved exactly in the notes: RUL1, finite double shuffle + duality);
# (2) the table of 307 monomials is folded exactly at x = -1 with these values; the coefficients of the 20 constants are compared with
#     the integer relation found numerically at 250 and 320 digits by the referee of the notes (the relation is not used in the fold).
import os, pickle, sys, time, sympy as sp, mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
WD = os.environ['G9WORK']
l2, pi, z3, z5, Li4h, Li5h = sp.symbols('l2 pi z3 z5 Li4h Li5h')
def Z23(c): return c*pi**2*z3/6
R = sp.Rational
atomval = {
 'H1': -l2, 'Li2': -pi**2/12, 'Li3': -3*z3/4, 'H011': z3/8, 'Li4': -7*pi**4/720,
 'H0011': -pi**4/48 + R(7, 4)*z3*l2 - R(1, 12)*pi**2*l2**2 + R(1, 12)*l2**4 + 2*Li4h,
 'H0111': -pi**4/90 + R(7, 8)*z3*l2 - R(1, 24)*pi**2*l2**2 + R(1, 24)*l2**4 + Li4h,
 'Li5': -R(15, 16)*z5,
 'H00011': -R(29, 32)*z5 + Z23(R(1, 2)),
 'H00101': R(41, 32)*z5 - Z23(R(5, 8)),
 'H00111': R(33, 32)*z5 + Z23(R(1, 2)) - R(7, 8)*z3*l2**2 + R(1, 18)*pi**2*l2**3 - R(1, 15)*l2**5 - 2*Li5h - 2*Li4h*l2,
 'H01011': -R(73, 64)*z5 - Z23(R(3, 2)) + R(7, 4)*z3*l2**2 - R(1, 9)*pi**2*l2**3 + R(2, 15)*l2**5 + 4*Li5h + 4*Li4h*l2,
 'H01111': z5 - R(7, 16)*z3*l2**2 + R(1, 36)*pi**2*l2**3 - R(1, 30)*l2**5 - Li5h - Li4h*l2,
}
# ---- (1) numerical check of the 13 values, x = -1 ----
mp.mp.dps = 110
NAMEWORD = {'H1': '1', 'Li2': '01', 'Li3': '001', 'H011': '011', 'Li4': '0001', 'H0011': '0011', 'H0111': '0111', 'Li5': '00001',
            'H00011': '00011', 'H00101': '00101', 'H00111': '00111', 'H01011': '01011', 'H01111': '01111'}
def mpl_indices(word):
    s, run = [], 0
    for ch in word:
        if ch == '0': run += 1
        else: s.append(run+1); run = 0
    return tuple(s)
def terms_nested(ss, N):                       # b_n = (1/n^{s1}) * sum_{n>n2>...} ..., n = 1..N
    arr = [mp.mpf(1)]*(N+2)
    for s in reversed(ss[1:]):
        new = [mp.mpf(0)]*(N+2); acc = mp.mpf(0)
        for j in range(1, N+1):
            acc += arr[j]/mp.mpf(j)**s
            new[j+1] = acc
        arr = new
    return [None] + [arr[n]/mp.mpf(n)**ss[0] for n in range(1, N+1)]
def cvz(a, N):                                 # sum_{k>=0} (-1)^k a_k  (Cohen, Villegas, Zagier), a_k = a[k], k < N
    d = (3 + mp.sqrt(8))**N; d = (d + 1/d)/2; b = mp.mpf(-1); c = -d; s = mp.mpf(0)
    for k in range(N):
        c = b - c; s += c*a[k]
        b = (k+N)*(k-N)*b/((k+mp.mpf(1)/2)*(k+1))
    return s/d
N = 170
cv = {'l2': mp.log(2), 'pi': mp.pi, 'z3': mp.zeta(3), 'z5': mp.zeta(5), 'Li4h': mp.polylog(4, mp.mpf(1)/2), 'Li5h': mp.polylog(5, mp.mpf(1)/2)}
worst = mp.mpf(0)
for nm, w in NAMEWORD.items():
    b = terms_nested(mpl_indices(w), N+2)
    ser = -cvz([b[k+1] for k in range(N)], N)            # sum_{n>=1} (-1)^n b_n
    ex = sp.lambdify(list(cv), atomval[nm], 'mpmath')(*[cv[k] for k in cv])
    worst = max(worst, abs(ser - ex))
print('13 harmonic polylogarithms at x = -1: worst |series (CVZ, %d terms) - closed form in the constants| = %s' % (N, mp.nstr(worst, 3)))
# ---- (2) exact fold ----
tab = pickle.load(open(f'{WD}/coef_acc_raw.pkl', 'rb'))
A0, A1, B, kap, L2, q0, V02, lA0, lA1, lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
z2s, z3s, z5s = sp.symbols('z2 z3 z5')
rd = {A0: R(1, 2), A1: R(3, 4), B: R(1, 4), kap: 1, L2: 5, q0: 1, V02: 1, lA0: -l2, lA1: -2*l2, lB: -2*l2}
t0 = time.time(); fold = sp.Integer(0)
for (at, xe, ze), co in tab.items():
    mv = sp.Integer((-1)**xe)
    for nm, e in at: mv *= atomval[nm]**e
    mv *= (pi**2/6)**ze[0]*z3**ze[1]*z5**ze[2]
    fold += co.subs(rd)*mv
fold = sp.expand(120*fold)
X = sp.expand(fold/sp.sqrt(5))
names = ['1', 'l2', 'pi2', 'l2^2', 'z3', 'pi2 l2', 'l2^3', 'pi4', 'z3 l2', 'pi2 l2^2', 'l2^4', 'Li4h', 'z5', 'z2z3', 'pi4 l2', 'z3 l2^2', 'pi2 l2^3', 'l2^5', 'Li5h', 'Li4h l2']
mons = [1, l2, pi**2, l2**2, z3, pi**2*l2, l2**3, pi**4, z3*l2, pi**2*l2**2, l2**4, Li4h, z5, pi**2*z3/6, pi**4*l2, z3*l2**2, pi**2*l2**3, l2**5, Li5h, Li4h*l2]
P = sp.Poly(X, l2, pi, z3, z5, Li4h, Li5h)
fitted = {}
for nm, m in zip(names, mons):
    mm = sp.Poly(m, l2, pi, z3, z5, Li4h, Li5h)
    (mon, _), = mm.terms()
    c = P.coeff_monomial(mon)*(mm.terms()[0][1]**-1)
    fitted[nm] = sp.nsimplify(c)
refrel = {'1': 0, 'l2': 96, 'pi2': R(-32, 5), 'l2^2': R(-576, 5), 'z3': R(3267, 50), 'pi2 l2': R(1287, 50), 'l2^3': R(222, 25), 'pi4': R(-6489, 8000),
          'z3 l2': R(-558, 5), 'pi2 l2^2': R(-4689, 200), 'l2^4': R(-39, 10), 'Li4h': R(-1242, 25), 'z5': R(-20187, 2000), 'z2z3': R(-858, 125),
          'pi4 l2': R(289, 500), 'z3 l2^2': R(6651, 200), 'pi2 l2^3': R(3309, 500), 'l2^5': R(2257, 500), 'Li5h': R(396, 125), 'Li4h l2': R(594, 25)}
print('exact fold (%.0f s): coefficients of E_e/sqrt5 on the 20 constants:' % (time.time()-t0))
print({k: str(v) for k, v in fitted.items()})
print('all 20 coefficients equal the integer relation of the referee (weighted convention):', all(sp.simplify(fitted[k] - refrel[k]) == 0 for k in names))
print('[zeta5] = %s sqrt5 ; [zeta2 zeta3] = %s sqrt5 ; remainder outside the span of the 20 constants: %s' % (fitted['z5'], fitted['z2z3'], sp.simplify(X - sum(fitted[k]*m for k, m in zip(names, mons)))))
# numerical value of the edge
val = sp.lambdify(list(cv), fold, 'mpmath')(*[cv[k] for k in cv])
print('E_e (r=2 edge, weighted convention) =', mp.nstr(val, 40))
pickle.dump(str(fold), open(f'{WD}/r2_fold.pkl', 'wb'))

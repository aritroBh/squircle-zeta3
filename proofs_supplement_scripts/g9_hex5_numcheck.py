# g9_hex5_numcheck.py -- fresh numerical check of the closed form of Theorem 26 (own code).
# The table of 307 monomials (pickle in $G9WORK, rebuilt by g9_hex5_chain.sh) is evaluated at 60 digits for hexagons
# whose three ratios r satisfy 0 < r < 2 (|x| < 1, so that the 13 harmonic polylogarithms are summed from their nested series),
#   P^(5)(1) = 2 * sum_{3 edges} 120 * sum_monomials coef * x^k * atoms * zeta's.
# The independent value is the finite-difference value of g9_direct5.py (perimeter of K_p from the definition).
import os, sys, pickle, time, mpmath as mp, sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g9_direct5 as dr
mp.mp.dps = 60
WD = os.environ['G9WORK']
tab = pickle.load(open(f'{WD}/coef_acc_raw.pkl', 'rb'))
A0, A1, B, kap, L2, q0, V02, lA0, lA1, lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
syms = (A0, A1, B, kap, L2, q0, V02, lA0, lA1, lB)
t0 = time.time()
fun = {k: sp.lambdify(syms, c, 'mpmath') for k, c in tab.items()}
print('lambdified %d coefficients (%.0f s)' % (len(fun), time.time()-t0))
NAMEWORD = {'H1': '1', 'Li2': '01', 'Li3': '001', 'H011': '011', 'Li4': '0001', 'H0011': '0011', 'H0111': '0111', 'Li5': '00001',
            'H00011': '00011', 'H00101': '00101', 'H00111': '00111', 'H01011': '01011', 'H01111': '01111'}
def mpl_indices(word):                         # 0^{s1-1} 1 0^{s2-1} 1 ...  ->  (s1, s2, ...)
    s, run = [], 0
    for ch in word:
        if ch == '0': run += 1
        else: s.append(run+1); run = 0
    return tuple(s)
def hpl(name, x, N):
    ss = mpl_indices(NAMEWORD[name])
    arr = [mp.mpf(1)]*(N+2)                     # tail = empty: nested sum 1
    for s in reversed(ss[1:]):                  # build the nested sums of the tail, innermost last
        new = [mp.mpf(0)]*(N+2); acc = mp.mpf(0)
        for j in range(1, N+1):
            acc += arr[j]/mp.mpf(j)**s
            new[j+1] = acc
        arr = new
    tot, xp = mp.mpf(0), mp.mpf(1)
    for n in range(1, N+1):
        xp *= x
        tot += xp*arr[n]/mp.mpf(n)**ss[0]
    return tot
# self-test of the series evaluator against integrals:  H011 = (1/2) int_0^x log^2(1-s)/s ds,  H1 = -log(1-x),  Li3 = polylog
xt = mp.mpf('0.3')
print('self-test H1:', mp.nstr(hpl('H1', xt, 300) + mp.log(1-xt), 3), ' Li3:', mp.nstr(hpl('Li3', xt, 300) - mp.polylog(3, xt), 3),
      ' H011:', mp.nstr(hpl('H011', xt, 300) - mp.quad(lambda s: mp.log(1-s)**2/s/2, [0, xt]), 3),
      ' H00101 (Li_{3,2}):', mp.nstr(hpl('H00101', xt, 300) - sum(xt**n/mp.mpf(n)**3*sum(1/mp.mpf(m)**2 for m in range(1, n)) for n in range(1, 300)), 3))
zv = (mp.zeta(2), mp.zeta(3), mp.zeta(5))
def edge_value(data, x):
    N = int(mp.ceil(65/(-mp.log10(abs(x))))) + 10
    H = {nm: hpl(nm, x, N) for nm in NAMEWORD}
    args = [data[s] for s in ('A0', 'A1', 'B', 'kap', 'L2', 'q0', 'V02', 'lA0', 'lA1', 'lB')]
    tot = mp.mpf(0)
    for k, f in fun.items():
        v = f(*args)*x**k[1]
        for nm, e in k[0]: v *= H[nm]**e
        for z, e in zip(zv, k[2]): v *= z**e
        tot += v
    return 120*tot
def hexagon_edges(C, ws):
    us = dr.us_from_corners(C)
    pts = list(C) + [(-a, -b) for (a, b) in C]
    res = []
    for k in range(3):
        V0, V1 = pts[k], pts[k+1]
        j0, j1 = k, (k+1) % 3; jf = 3 - j0 - j1
        a = lambda j, V: abs(us[j][0]*V[0] + us[j][1]*V[1])
        Av0, Av1, Bv = a(j0, V1), a(j1, V0), a(jf, V0)
        r = a(jf, V1)/Bv
        assert abs(Av1 + Bv - 1) < mp.mpf(10)**-50 and abs(Av0 + Bv*r - 1) < mp.mpf(10)**-50
        D = (V1[0]-V0[0], V1[1]-V0[1])
        data = dict(A0=Av0, A1=Av1, B=Bv, kap=r-1, L2=D[0]**2 + D[1]**2, q0=V0[0]*D[0] + V0[1]*D[1], V02=V0[0]**2 + V0[1]**2,
                    lA0=mp.log(Av0/ws[j0]), lA1=mp.log(Av1/ws[j1]), lB=mp.log(Bv/ws[jf]))
        res.append((data, 1 - r))
    return res
R = mp.mpf
cases = {
 'A': ([(R(1), R(1)), (R(-27)/10, R(33)/10), (R(-2), R(1))], [5, 1, 2], '-1474.6784800498412080730263501189100249635'),
 'C': ([(R(2), -R(1)/3), (R(26)/7, R(3)/14), (R(1), R(1))], [R(3)/2, 4, 1], '720.74638522850829292847185026162593439825'),
}
which = sys.argv[1]
C, ws, ref = cases[which]
edges = hexagon_edges(C, ws)
print('hexagon', which, ': ratios r =', [mp.nstr(1-x, 8) for d, x in edges])
tot = 2*sum(edge_value(d, x) for d, x in edges)
print('P5(1) from the closed form  =', mp.nstr(tot, 40))
us = dr.us_from_corners(C)
T, P1 = dr.P5(us, ws, mp.mpf('0.002'), 7)
print('P5(1) finite differences     =', mp.nstr(T[-1][-1], 40), ' (dps %d, h0 0.002, 7 Richardson levels)' % mp.mp.dps)
print('reference of the notes (referee, 110 digits) =', ref[:42])
print('|closed form - finite differences| =', mp.nstr(abs(tot - T[-1][-1]), 3), '   |closed form - reference| =', mp.nstr(abs(tot - mp.mpf(ref)), 3))
print('done (%.0f s)' % (time.time()-t0))

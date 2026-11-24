# g9_hex5_r1limit.py -- the table of Theorem 26 at r = 1 (kappa -> 0), against the exact local value of the same edge (own code).
# Hexagon (1,0),(2,2),(0,1), weights (1,3,1): its third edge has r = 1.  The table (pickle in $G9WORK) is evaluated at x = +-10^-12 with
# the other edge data fixed (A0 = 1 - B - B kappa), at 150 digits (the monomials have poles x^-6), and averaged over the two signs.
# The local value is E_e = 120 * int E_5 dt of the exact fully local engine of g9_local5.py (constant far term B, logarithms lA0, lA1, lB).
import os, sys, pickle, mpmath as mp, sympy as sp
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g9_local5 as loc          # prints its own checks when imported
mp.mp.dps = 150
WD = os.environ['G9WORK']
tab = pickle.load(open(os.environ.get('G9TABLE', f'{WD}/coef_acc_raw.pkl'), 'rb'))
A0, A1, B, kap, L2, q0, V02, lA0, lA1, lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
syms = (A0, A1, B, kap, L2, q0, V02, lA0, lA1, lB)
fun = {k: sp.lambdify(syms, c, 'mpmath') for k, c in tab.items()}
NAMEWORD = {'H1': '1', 'Li2': '01', 'Li3': '001', 'H011': '011', 'Li4': '0001', 'H0011': '0011', 'H0111': '0111', 'Li5': '00001',
            'H00011': '00011', 'H00101': '00101', 'H00111': '00111', 'H01011': '01011', 'H01111': '01111'}
def mpl_indices(word):
    s, run = [], 0
    for ch in word:
        if ch == '0': run += 1
        else: s.append(run+1); run = 0
    return tuple(s)
def hpl(name, x, N):
    ss = mpl_indices(NAMEWORD[name]); arr = [mp.mpf(1)]*(N+2)
    for s in reversed(ss[1:]):
        new = [mp.mpf(0)]*(N+2); acc = mp.mpf(0)
        for j in range(1, N+1):
            acc += arr[j]/mp.mpf(j)**s; new[j+1] = acc
        arr = new
    tot, xp = mp.mpf(0), mp.mpf(1)
    for n in range(1, N+1):
        xp *= x; tot += xp*arr[n]/mp.mpf(n)**ss[0]
    return tot
zv = (mp.zeta(2), mp.zeta(3), mp.zeta(5))
def table_value(d, x):
    H = {nm: hpl(nm, x, 20) for nm in NAMEWORD}
    args = [d['A0'], d['A1'], d['B'], -x, d['L2'], d['q0'], d['V02'], d['lA0'], d['lA1'], d['lB']]
    args[0] = 1 - d['B'] - d['B']*(-x)      # A0 = 1 - B - B kappa, kappa = -x
    tot = mp.mpf(0)
    for k, f in fun.items():
        v = f(*args)*x**k[1]
        for nm, e in k[0]: v *= H[nm]**e
        for z, e in zip(zv, k[2]): v *= z**e
        tot += v
    return 120*tot
# exact edge data of the third edge of the hexagon (1,0),(2,2),(0,1) (corners C0,C1,C2, then -C0): edge from C2 to -C0
R = sp.Rational
C = [(R(1), R(0)), (R(2), R(2)), (R(0), R(1))]; w = [1, 3, 1]
def grad(P, Q):
    det = P[0]*Q[1] - P[1]*Q[0]; return ((Q[1]-P[1])/det, (P[0]-Q[0])/det)
pts = C + [(-a, -b) for a, b in C]
g = [grad(pts[k], pts[(k+1) % 6]) for k in range(6)]
u = [((g[k][0]-g[k-1][0])/2, (g[k][1]-g[k-1][1])/2) for k in range(3)]
V0, V1 = pts[2], pts[3]
a = lambda j, V: abs(u[j][0]*V[0] + u[j][1]*V[1])
j0, j1, jf = 2, 0, 1
alpha, beta, Bv = a(j0, V1), a(j1, V0), a(jf, V0)
r = a(jf, V1)/Bv
D = (V1[0]-V0[0], V1[1]-V0[1]); L2v = D[0]**2 + D[1]**2; q0v = V0[0]*D[0] + V0[1]*D[1]; V02v = V0[0]**2 + V0[1]**2
print('third edge of the paper hexagon: A0, A1, B, r =', alpha, beta, Bv, r, ' L2, q0, V02 =', L2v, q0v, V02v)
assert r == 1 and alpha + Bv == 1 and beta + Bv == 1
logs = dict(lA0=mp.log(mp.mpf(alpha.p)/alpha.q/w[j0]), lA1=mp.log(mp.mpf(beta.p)/beta.q/w[j1]), lB=mp.log(mp.mpf(Bv.p)/Bv.q/w[jf]))
data = dict(A0=mp.mpf(alpha.p)/alpha.q, A1=mp.mpf(beta.p)/beta.q, B=mp.mpf(Bv.p)/Bv.q, L2=mp.mpf(L2v.p)/L2v.q, q0=mp.mpf(q0v.p)/q0v.q, V02=mp.mpf(V02v.p)/V02v.q, **logs)
x0 = mp.mpf(10)**-12
vp, vm = table_value(data, x0), table_value(data, -x0)
tab_lim = (vp + vm)/2
print('table at x = +10^-12 :', mp.nstr(vp, 30)); print('table at x = -10^-12 :', mp.nstr(vm, 30))
# exact local value of the same edge
u0 = q0v/L2v; mm = (V02v - q0v**2/L2v)/L2v
mu0, mu1, lbs = sp.symbols('mu0 mu1 lbs')
I = loc.edge_integrals(alpha, beta, [(Bv, lbs)], mu0, mu1, u0, mm)
Lv = mp.sqrt(data['L2'])
val = sp.lambdify((mu0, mu1, lbs, loc.z2, loc.z3, loc.z5), I[5], 'mpmath')(logs['lA0'], logs['lA1'], logs['lB'], mp.zeta(2), mp.zeta(3), mp.zeta(5))
local = 120*Lv*val
print('exact local value of the edge (g9_local5 engine):', mp.nstr(local, 30))
print('|table limit - local value| =', mp.nstr(abs(tab_lim - local), 3), '   (error of the average of the two signs is O(x^2) = 1e-24)')

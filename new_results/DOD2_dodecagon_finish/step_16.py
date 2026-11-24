
import mpmath as mp, pickle, time
import sympy as sp
from sympy import Rational, sqrt, symbols, expand, Function, simplify, Add, Mul
mp.mp.dps = 300
t0=time.time()
r3 = sqrt(3)
l2,l3,lu = symbols('l2 l3 lu')
Li2A = Function('Li2A'); Li3A = Function('Li3A')
Xs = {1: 1-2*r3/3, 2: 1-r3, 5: 1-r3/3, 6: 1-r3/2}
UP = 2+r3
import mpmath as mpx
def logfac(theta):
    th = sp.simplify(sp.expand(theta))
    if sp.simplify(th) == 1: return sp.Integer(0)
    sgn = 1
    if th < 0: sgn = -1; th = -th
    mpx.mp.dps = 100
    lv = mpx.mpf(sp.N(sp.log(th), 90)); L2 = mpx.log(2); L3 = mpx.log(3); LU = mpx.log(2-mpx.sqrt(3))
    rel = mpx.pslq([lv, L2, L3, LU], tol=mpx.mpf('1e-85'), maxcoeff=10**7, maxsteps=4000)
    assert rel is not None and rel[0] != 0, "no logfac for %s"%theta
    a = -Rational(rel[1], rel[0]); b = -Rational(rel[2], rel[0]); c = -Rational(rel[3], rel[0])
    cand = 2**a * 3**b * (2-r3)**c
    assert sp.simplify(sp.expand(th - cand)) == 0, "logfac verify fail %s"%theta
    return a*l2 + b*l3 + c*lu
LOG1 = {i: logfac(1-Xs[i]) for i in Xs}
def L1abs(kap,u):
    # returns sympy expr; kap,u exact in Q(sqrt3); branch on kap*u vs 1
    if sp.simplify(kap*u - 1) <= 0:
        return -Li2A(kap*u)
    else:
        return logfac(kap)*logfac(u) + logfac(u)**2/2 + Li2A(1/(kap*u))
def I2abs(kap,u):
    if sp.simplify(kap*u - 1) <= 0:
        return -logfac(u)*Li2A(kap*u) + Li3A(kap*u)
    else:
        return logfac(kap)*logfac(u)**2/2 + logfac(u)**3/3 + logfac(u)*Li2A(1/(kap*u)) + Li3A(1/(kap*u))
def Mcf(c, xi):
    kap = c/(c-1); lo, hi = 1-xi, sp.Integer(1)
    return ( logfac(abs(1-c))*(logfac(hi)**2-logfac(lo)**2)/2 + I2abs(kap,hi) - I2abs(kap,lo)
             - (logfac(hi)**3-logfac(lo)**3)/3 )
def QnonC(a,b,xi):
    ka, kb = a/(a-1), b/(b-1)
    A, B = logfac(abs(1-a)), logfac(abs(1-b))
    lo, hi = 1-xi, sp.Integer(1)
    def F(u):
        return ( A*B*logfac(u) + A*L1abs(kb,u) - A*logfac(u)**2/2
               + B*L1abs(ka,u) - B*logfac(u)**2/2
               - I2abs(ka,u) - I2abs(kb,u) + logfac(u)**3/3 )
    return F(hi)-F(lo)
def Swrem(ka,kb):
    return -logfac(abs(1-ka))*(Li2A(kb)-Li2A(kb/ka)) + Li3A(kb)
i,j,k = 2,5,6
xi = Xs[i]; a = Xs[j]/xi; b = Xs[k]/xi
ka, kb = a/(a-1), b/(b-1)
lo = 1-xi
da, db = ka*lo, kb*lo
bnd = -logfac(1-xi)**3/3
R = (bnd + Mcf(a,xi) + Mcf(b,xi) + QnonC(a,b,xi) + Swrem(ka,kb) - Swrem(da,db))/xi
R = expand(R)
print("R built, terms:", len(Add.make_args(R)), "t=%ds"%(time.time()-t0))
# collect atoms
atoms = R.atoms(Li2A) | R.atoms(Li3A)
print("polylog atoms:", sorted([str(s) for s in atoms]))
pickle.dump(R, open('R256_sym.pkl','wb'))
# numeric eval
mp.mp.dps = 620
vals = {l2: mp.log(2), l3: mp.log(3), lu: mp.log(2-mp.sqrt(3))}
def ev(e):
    e = e.subs(vals)
    def r2(m):
        f = m.func
        if f == Li2A: return mp.re(mp.polylog(2, mp.mpc(sp.N(m.args[0], 700))))
        if f == Li3A: return mp.re(mp.polylog(3, mp.mpc(sp.N(m.args[0], 700))))
        return m
    return sp.N(e.replace(Li2A, lambda z: mp.re(mp.polylog(2, mp.mpc(sp.N(z,700))))).replace(Li3A, lambda z: mp.re(mp.polylog(3, mp.mpc(sp.N(z,700))))), 620)
import sympy
Rv = ev(R)
D1 = '../DOD1_dodecagon/'
atoms620 = pickle.load(open(D1+'atoms620.pkl','rb'))
g = pickle.load(open('gdswaps620.pkl','rb'))
R_quad = atoms620['T3a'][(i,j,k)] - (g['SwG'][(i,j,k)]-g['SwD'][(i,j,k)])/atoms620['X'][i]
print("R_sym =", mp.nstr(mp.mpf(Rv), 20))
print("R_quad=", mp.nstr(R_quad, 20))
print("|diff| =", mp.nstr(abs(mp.mpf(Rv)-R_quad),4))
print("elapsed", time.time()-t0)

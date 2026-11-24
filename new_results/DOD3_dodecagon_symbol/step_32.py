import sys, os, pickle
sys.path.insert(0, os.getcwd())
import mpmath as mp
mp.mp.dps = 200
import sympy as sp
from dod3_R import *

# numeric piecewise definitions (step_14 style, with Cre by quadrature)
Li2R = lambda z: mp.re(mp.polylog(2, mp.mpc(z)))
Li3R = lambda z: mp.re(mp.polylog(3, mp.mpc(z)))
def L1absN(kap,u):
    if kap*u <= 1: return -Li2R(kap*u)
    else: return mp.log(kap)*mp.log(u) + mp.log(u)**2/2 + Li2R(1/(kap*u))
def I2absN(kap,u):
    if kap*u <= 1: return -mp.log(u)*Li2R(kap*u) + Li3R(kap*u)
    else: return mp.log(kap)*mp.log(u)**2/2 + mp.log(u)**3/3 + mp.log(u)*Li2R(1/(kap*u)) + Li3R(1/(kap*u))
def CreN(al,be,u):
    pts = [mp.mpf(0)]
    for kap in (al,be):
        if kap > 0:
            s0 = 1/kap
            if 0 < s0 < u: pts.append(s0)
    pts.append(u); pts = sorted(set(pts))
    tot = mp.mpf(0)
    for a_,b_ in zip(pts[:-1], pts[1:]):
        tot += mp.quad(lambda s: mp.log(abs(1-al*s))*mp.log(abs(1-be*s))/s, [a_,b_])
    return tot
def McN(c, xi):
    kap = c/(c-1); lo, hi = 1-xi, mp.mpf(1)
    return ( mp.log(abs(1-c))*(mp.log(hi)**2-mp.log(lo)**2)/2 + I2absN(kap,hi) - I2absN(kap,lo)
             - (mp.log(hi)**3-mp.log(lo)**3)/3 )
def QabN(a,b,xi):
    ka, kb = a/(a-1), b/(b-1)
    A, B = mp.log(abs(1-a)), mp.log(abs(1-b))
    lo, hi = 1-xi, mp.mpf(1)
    def F(u):
        return ( A*B*mp.log(u) + A*L1absN(kb,u) - A*mp.log(u)**2/2
               + B*L1absN(ka,u) - B*mp.log(u)**2/2 + CreN(ka,kb,u)
               - I2absN(ka,u) - I2absN(kb,u) + mp.log(u)**3/3 )
    return F(hi)-F(lo)
def SwpRN(a,b):
    return CreN(a,b,mp.mpf(1)) + mp.log(abs(1-a))*(Li2R(b)-Li2R(b/a)) - Li3R(b)

# symbolic evaluator
vals200 = {l2: mp.log(2), l3: mp.log(3), lu: mp.log(2-mp.sqrt(3))}
def ev(e):
    e = e.subs(vals200)
    e = e.replace(Li2A, lambda z: mp.re(mp.polylog(2, mp.mpc(sp.N(z,300))))).replace(Li3A, lambda z: mp.re(mp.polylog(3, mp.mpc(sp.N(z,300)))))
    return mp.mpf(sp.N(e, 200))

for (i,j,k) in [(6,1,5),(1,2,5),(2,5,6)]:
    xi = Xs[i]; a = pdiv(Xs[j],xi); b = pdiv(Xs[k],xi)
    an, bn, xin = pfloat(a), pfloat(b), pfloat(xi)
    kan, kbn = an/(an-1), bn/(bn-1)
    lo = 1-xin
    # numeric J pieces
    Jn = -mp.log(1-xin)**3/3 + McN(an,xin) + McN(bn,xin) + QabN(an,bn,xin)
    # my symbolic pieces
    Mc_a_s = Mcf(a,xi); Mc_b_s = Mcf(b,xi); Q_s = QnonC(a,b,xi)
    bnd_s = -logfac_pair(psub(ONE,xi))**3/3
    print("="*40, (i,j,k))
    print("  bnd:   sym %s  num %s"%(mp.nstr(ev(bnd_s),12), mp.nstr(-mp.log(1-xin)**3/3,12)))
    print("  Mc(a): sym %s  num %s"%(mp.nstr(ev(Mc_a_s),12), mp.nstr(McN(an,xin),12)))
    print("  Mc(b): sym %s  num %s"%(mp.nstr(ev(Mc_b_s),12), mp.nstr(McN(bn,xin),12)))
    # QnonC symbolic vs QabN - CreN + (C-C) swap... QnonC = Qab - Cre; numeric Q without Cre:
    Qnc_n = QabN(an,bn,xin) - (CreN(kan,kbn,mp.mpf(1)) - CreN(kan,kbn,lo))
    print("  QnonC: sym %s  num %s"%(mp.nstr(ev(Q_s),12), mp.nstr(Qnc_n,12)))

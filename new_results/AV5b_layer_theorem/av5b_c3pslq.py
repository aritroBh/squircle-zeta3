
import mpmath as mp
import json
mp.mp.dps = 80
from scipy.interpolate import CubicSpline
import numpy as np

def load_prof(mm):
    return [(mp.mpf(s), mp.mpf(a), mp.mpf(b), mp.mpf(c)) for s,a,b,c in json.load(open(f"ckpt_prof_m{mm}.json"))]

def c3_lay_accurate(mm):
    prof = load_prof(mm)
    ss = [p[0] for p in prof]
    H = mp.sin(mp.pi/(2*mm)); alpha = mp.tan(mp.pi/(2*mm))
    a1 = prof[-1][1]; a2 = prof[-1][2]
    # extend tails analytically: g1 ~ A e^{-r|s|} fits at the last two points for quad tails
    g1v = [float(p[1]-a1) for p in prof]; g2v = [float(p[2]-a2) for p in prof]
    sf = [float(s) for s in ss]
    sp1 = CubicSpline(sf, g1v, bc_type='natural')
    sp2 = CubicSpline(sf, g2v, bc_type='natural')
    g1 = lambda s: mp.mpf(str(sp1(float(s))))
    g2 = lambda s: mp.mpf(str(sp2(float(s))))
    I_g2 = mp.quad(g2, [-20, 0, 20])
    I_g1sq = mp.quad(lambda s: g1(s)**2/2, [-20, 0, 20])
    I_sg1 = mp.quad(lambda s: (alpha*abs(s)+a1)*g1(s), [-20, 0, 20])
    return H*(I_g2 + I_g1sq + I_sg1), (I_g2, I_g1sq, I_sg1), a1, a2

for mm, scale, basis_extra in [(3, mp.sqrt(3), []), (4, mp.sqrt(2+mp.sqrt(2)), [])]:
    c3, pieces, a1, a2 = c3_lay_accurate(mm)
    print(f"m={mm}: c3_lay = {mp.nstr(c3, 40)}")
    print(f"   pieces: I_g2={mp.nstr(pieces[0],18)} I_g1sq={mp.nstr(pieces[1],18)} I_sg1={mp.nstr(pieces[2],18)}")
    z3 = mp.zeta(3); ln2 = mp.log(2); pi2 = mp.pi**2
    basis = [z3, z3*ln2, pi2*ln2, ln2**3, pi2, ln2**2, ln2, mp.mpf(1)]
    M = [c3/scale] + basis
    rel = mp.pslq(M, tol=mp.mpf('1e-25'), maxcoeff=10**12, maxsteps=10000)
    print(f"   PSLQ c3/scale:", rel)
    if rel:
        print("    resid:", mp.nstr(mp.fsum(mp.mpf(r)*x for r,x in zip(rel,M)),6))
        # zeta(3) coefficient:
        print("    => zeta(3) coeff in c3_lay =", -rel[1]/rel[0], "* scale")

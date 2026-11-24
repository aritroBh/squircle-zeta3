import sys, os
workdir = '<project>/research/08_theory/TOW2R_review'
sys.path.insert(0, workdir); sys.path.insert(0, os.path.join(workdir, 'tow2_copy'))
import tow2r_lib as RL
import buildintegrand as BI
import genintegrand as GI
import fastred
from fastred import QSI
import mpmath as mp
mp.mp.dps = 40
def qsi_to_mp(z):
    return mp.mpf(z.p0.numerator)/z.p0.denominator + mp.mpf(z.p1.numerator)/z.p1.denominator*mp.sqrt(2) \
           + 1j*(mp.mpf(z.p2.numerator)/z.p2.denominator) + 1j*mp.mpf(z.p3.numerator)/z.p3.denominator*mp.sqrt(2)
def peval_mp(P, v):
    out = mp.mpc(0)
    for c in reversed(P): out = out*v + qsi_to_mp(c)
    return out

for (p,q) in [(4,0),(3,1),(2,2)]:
    num, den = RL.Phi_num_den(p, q)
    S = GI.S_poly(p+q, p)
    PF = BI.sqrt_part_PF(p+q, S)
    for vx, vm in [(QSI(1,0,1,0), 1+1j), (QSI(2,0,0,1), 2+mp.sqrt(2)*1j)]:
        mine2 = peval_mp(num, vm)/peval_mp(den, vm)
        ref2 = qsi_to_mp(fastred.pf_eval_at(PF, vx))
        print((p,q), 'ratio mine/ref =', mp.nstr(mine2/ref2, 8) if ref2 else 'ref=0')

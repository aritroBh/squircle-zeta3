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

num, den = RL.Phi_num_den(4, 0)
v = mp.mpc('0.37', '0.11')
mine = peval_mp(num, v)/peval_mp(den, v)
S = GI.S_poly(4, 4)
PF = BI.sqrt_part_PF(4, S)
ref = qsi_to_mp(fastred.pf_eval_at(PF, QSI(0)) ) # need x as QSI: v=0.37+0.11i not in Q(sqrt2,i)... use exact point
# use an exact QSI point instead
vx = QSI(1, 0, 1, 0)  # 1 + i
mine2 = peval_mp(num, 1+1j)/peval_mp(den, 1+1j)
ref2 = qsi_to_mp(fastred.pf_eval_at(PF, vx))
print('mine:', mine2)
print('ref :', ref2)
print('ratio:', mine2/ref2 if ref2 else 'ref=0')

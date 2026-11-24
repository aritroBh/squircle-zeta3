import sys, os
workdir = '<project>/research/08_theory/TOW2R_review'
sys.path.insert(0, workdir); sys.path.insert(0, os.path.join(workdir, 'tow2_copy'))
import tow2r_lib as RL
import buildintegrand as BI
import genintegrand as GI
import mpmath as mp
mp.mp.dps = 40

# my Phi for (4,0) as num/den, evaluate at v
num, den = RL.Phi_num_den(4, 0)
def qsi_to_mp(z):
    return mp.mpf(z.p0.numerator)/z.p0.denominator + mp.mpf(z.p1.numerator)/z.p1.denominator*mp.sqrt(2) \
           + 1j*(mp.mpf(z.p2.numerator)/z.p2.denominator) + 1j*mp.mpf(z.p3.numerator)/z.p3.denominator*mp.sqrt(2)
def peval_mp(P, v):
    out = mp.mpc(0)
    for c in reversed(P): out = out*v + qsi_to_mp(c)
    return out
v = mp.mpc('0.37', '0.11')
mine = peval_mp(num, v)/peval_mp(den, v)

# buildintegrand sqrt_part_PF for S_poly(4,4) (m=p=4): unfolded piece
S = GI.S_poly(4, 4)
PF = BI.sqrt_part_PF(4, S)
print(type(PF), len(PF))
# how to evaluate PF? look at fastred
import inspect, fastred
src = inspect.getsource(fastred)
ix = src.find('def pf_eval')
print(src[ix:ix+400] if ix>=0 else 'no pf_eval; PF structure:')
ix = src.find('def pf_from_poly')
print(src[ix:ix+700])

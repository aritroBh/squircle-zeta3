import sys, os
workdir = '<project>/research/08_theory/TOW2R_review'
sys.path.insert(0, workdir); sys.path.insert(0, os.path.join(workdir, 'tow2_copy'))
import tow2r_lib as RL
import buildintegrand as BI
import genintegrand as GI
from fastred import QSI, ZERO, ONE

S = GI.S_poly(4, 4)
print('S coeffs:', S)
# direct check: Stilde(v) should equal S(t(v))*(4v)^D at a sample point
from fractions import Fraction as Fr
def qsi_peval(P, x):
    tot = ZERO
    for c in reversed(P): tot = tot*x + c
    return tot
vx = QSI(2)
tv = (vx*vx + 2*vx - 1)/(4*vx)
St_val = sum(QSI(c)*tv**k for k,c in enumerate(S)) * (4*vx)**(len(S)-1)
mine, _ = RL.Stilde(S); ref, _ = BI.Stilde(S)
print('S(t(v))*(4v)^D at v=2:', St_val)
print('mine at v=2          :', qsi_peval(mine, vx))
print('ref  at v=2          :', qsi_peval(ref, vx))

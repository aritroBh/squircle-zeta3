import sys, time; sys.dont_write_bytecode=True
sys.path.insert(0,'.')
from fractions import Fraction as F
import phi
t0=time.time()
for W in (11,):
    t=time.time(); phi.E._TABLES[W]=phi.M.DSTableFlint(W); print('table',W,'built %.1fs'%(time.time()-t),flush=True)
for a in [(1,1,1,2),(3,1,2,5)]:
    k=phi.kin_from_a(a); e0,e1=phi.e_mats(**k)
    out=phi.phi_first_column(e0,e1,11)
    eng=phi.E.reduced_coeffs(a,11)
    print(a,'W=11 match engine:', out[11][0]==eng[11],' zeta(3,3,5) coeff:',out[11][0][-1],' zeta(3,5)zeta(3) coeff:',out[11][0][-2],flush=True)
print('total %.1fs'%(time.time()-t0))

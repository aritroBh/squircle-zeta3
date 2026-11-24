import sys, time; sys.dont_write_bytecode=True
sys.path.insert(0,'.')
from fractions import Fraction as F
import phi
t0=time.time()
for W in (9,10):
    t=time.time(); phi.E._TABLES[W]=phi.M.DSTableFlint(W); print('table',W,'built %.1fs'%(time.time()-t),flush=True)
a=(1,1,1,2)
k=phi.kin_from_a(a); e0,e1=phi.e_mats(**k)
t=time.time()
out=phi.phi_first_column(e0,e1,10)
print('phi to weight 10 %.1fs'%(time.time()-t),flush=True)
eng=phi.E.reduced_coeffs(a,10)
for W in (9,10):
    print('W=%d match engine:'%W, out[W][0]==eng[W], ' basis:',phi.E.table(W).monomials)
print('W=10 zeta(3,7) coefficient (last basis elt) from associator:',out[10][0][-1],' (task says 29515/56)')
print('total %.1fs'%(time.time()-t0))

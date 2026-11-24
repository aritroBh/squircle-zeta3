import sys, time; sys.dont_write_bytecode=True
sys.path.insert(0,'.')
from fractions import Fraction as F
import phi
pts=[(1,1,1,2),(1,2,2,1),(3,1,2,5),(2,1,1,3)]
Wmax=8
for a in pts:
    t=time.time()
    k=phi.kin_from_a(a)
    e0,e1=phi.e_mats(**k)
    out=phi.phi_first_column(e0,e1,Wmax)
    eng=phi.E.reduced_coeffs(a,Wmax)
    print('a=',a,'kin',{x:str(y) for x,y in k.items()})
    for W in range(2,Wmax+1):
        mine=out[W][0]     # F^{(23)} = Phi_11 * F2, F2=1 at s23=0
        ok = (mine==eng[W])
        print('  W=%d basis=%s match engine: %s'%(W,[m for m in phi.E.table(W).monomials],ok))
        if W==8:
            print('   Phi_11 at W=8 (basis zeta8, z3z5, z2z3^2, zeta(3,5)) =',[str(x) for x in mine])
            print('   engine                                          =',[str(x) for x in eng[W]])
    print('  time %.1fs'%(time.time()-t))

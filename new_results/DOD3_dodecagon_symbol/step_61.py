import sys, os
sys.path.insert(0, os.getcwd())
import mpmath as mp
mp.mp.dps = 620
l2,l3,lu = mp.log(2),mp.log(3),mp.log(2+mp.sqrt(3)); z2,z3c = mp.pi**2/6, mp.zeta(3)
L3c = (mp.zeta(3, mp.mpf(1)/12) - mp.zeta(3, mp.mpf(5)/12) - mp.zeta(3, mp.mpf(7)/12) + mp.zeta(3, mp.mpf(11)/12)) / 12**3
L2c = (mp.zeta(2, mp.mpf(1)/12) - mp.zeta(2, mp.mpf(5)/12) - mp.zeta(2, mp.mpf(7)/12) + mp.zeta(2, mp.mpf(11)/12)) / 12**2
Li3_2 = mp.re(mp.polylog(3, mp.mpc(2)))
print("Re Li3(2) =", mp.nstr(Li3_2, 25))
# known closed form? test Li3(2) vs [z3, z2*l2]
rel = mp.pslq([Li3_2, z3c, z2*l2, l2**3], tol=mp.mpf('1e-590'), maxcoeff=10**10, maxsteps=10000)
print("ReLi3(2) vs [z3, z2 l2, l2^3]:", rel)
rem = mp.mpf('-32.799586563322169122')
basis = [z3c, L3c, z2*l2, z2*l3, z2*lu, L2c*l2, L2c*l3, L2c*lu, Li3_2]
names = ['z3','L3chi','z2l2','z2l3','z2lu','L2l2','L2l3','L2lu','ReLi3(2)']
rel2 = mp.pslq([rem]+basis, tol=mp.mpf('1e-580'), maxcoeff=10**12, maxsteps=50000)
print("PSLQ rem:", rel2)
if rel2:
    resid = abs(rel2[0]*rem + sum(rel2[i+1]*basis[i] for i in range(len(basis))))
    print("resid:", mp.nstr(resid,3))
    for n,r in zip(['rem']+names, rel2):
        if r: print("   %-10s %s"%(n,r))

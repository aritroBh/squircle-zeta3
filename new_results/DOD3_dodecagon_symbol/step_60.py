import sys, os
sys.path.insert(0, os.getcwd())
import mpmath as mp
mp.mp.dps = 620
l2,l3,lu = mp.log(2),mp.log(3),mp.log(2+mp.sqrt(3)); z2,z3c = mp.pi**2/6, mp.zeta(3)
L3c = (mp.zeta(3, mp.mpf(1)/12) - mp.zeta(3, mp.mpf(5)/12) - mp.zeta(3, mp.mpf(7)/12) + mp.zeta(3, mp.mpf(11)/12)) / 12**3
L2c = (mp.zeta(2, mp.mpf(1)/12) - mp.zeta(2, mp.mpf(5)/12) - mp.zeta(2, mp.mpf(7)/12) + mp.zeta(2, mp.mpf(11)/12)) / 12**2
print("L(2,chi_12) =", mp.nstr(L2c, 40))
# series check
mp.mp.dps = 100
chi = {1:1, 5:-1, 7:-1, 11:1}
s2 = mp.mpf(0)
for n in range(1, 2000001):
    c = chi.get(n%12, 0)
    if c: s2 += c/mp.power(n,2)
print("series check |diff|:", mp.nstr(abs(s2-L2c),3))
mp.mp.dps = 620
rem = mp.mpf('-32.799586563322169122')
basis = [z3c, L3c, z2*l2, z2*l3, z2*lu, L2c*l2, L2c*l3, L2c*lu]
names = ['z3','L3chi','z2l2','z2l3','z2lu','L2chi*l2','L2chi*l3','L2chi*lu']
rel = mp.pslq([rem]+basis, tol=mp.mpf('1e-590'), maxcoeff=10**12, maxsteps=30000)
print("PSLQ [rem, z3, L3chi, z2*logs, L2chi*logs]:", rel)
if rel:
    resid = abs(rel[0]*rem + sum(rel[i+1]*basis[i] for i in range(len(basis))))
    print("resid:", mp.nstr(resid,3))
    for n,r in zip(['rem']+names, rel):
        if r: print("   %-10s %s"%(n,r))

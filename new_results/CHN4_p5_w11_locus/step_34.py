import pickle
from fractions import Fraction as F
from flint import fmpq_mat, fmpq
d=pickle.load(open('interp_points.pkl','rb')); sel=[tuple(x) for x in d['sel']]; orbits=d['orbits']; mons7=d['mons7']
ev=pickle.load(open('f11_evals.pkl','rb'))
orbit_of={}
for i,o in enumerate(orbits):
    for e in o: orbit_of[e]=i
NORB=len(orbits)
def orbit_vals_int(a):
    vals=[0]*NORB
    for e in mons7:
        pw=1
        for x,ei in zip(a,e):
            if ei: pw*=x**ei
        vals[orbit_of[e]]+=pw
    return vals
def Lval(a):
    a1,a2,a3,a4,a5=a
    return a2*a3*a4*(a1+a2+a3+a4+a5)

# build system
fitpts=sel
Mdata=[]; rhs=[]
for a in fitpts:
    lb,v=ev[a]; assert lb=='fit'
    Mdata += [fmpq(v) for v in orbit_vals_int(a)]
    rhs.append(fmpq(v.numerator, v.denominator)/Lval(a))
Mm=fmpq_mat(NORB,NORB,Mdata)
bm=fmpq_mat(NORB,1,rhs)
import time; t0=time.time()
sol=Mm.solve(bm)
print(f"solve took {time.time()-t0:.1f}s")
coefs=[F(int(sol[i,0].p), int(sol[i,0].q)) for i in range(NORB)]
# denominator analysis
from math import gcd
dens=[c.denominator for c in coefs]
lcm=1
for x in dens: lcm=lcm*x//gcd(lcm,x)
print("number of nonzero orbit coefs:", sum(1 for c in coefs if c!=0), "of", NORB)
print("lcm of denominators:", lcm)
print("coefs after clearing lcm are integers:", all((c*lcm).denominator==1 for c in coefs))
pickle.dump({'coefs':coefs,'orbits':orbits,'mons7':mons7,'lcm':lcm}, open('Q11_fit.pkl','wb'))
print("saved Q11_fit.pkl")

import os, time
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
import mpmath
mpmath.mp.dps = 120
from mpmath import mpf, sqrt, log, pi, zeta, quad
import sys; sys.path.insert(0,'.')
from lead2_lib import E_coeffs, m_e, L1, L2
from math import comb
from mpmath import exp as mexp
E = E_coeffs(20)
# (8): A_{p,q} = M(p,q;1/(1-u)) == Bose form
worst=mpf(0)
for (p,q) in [(1,0),(2,1),(3,2),(5,3),(4,4)]:
    M = quad(lambda u: L1(u)**p*L2(u)**q/(1-u), [-1,0,1])
    B = (-1)**p*quad(lambda r: r**p*log(1-mexp(-r))**q/(mexp(r)-1), [0,1,mpmath.inf])
    worst=max(worst,abs(M-B))
print("(8) M(p,q;1/(1-u)) == Bose form, max |diff| (p,q) in {(1,0),(2,1),(3,2),(5,3),(4,4)}:", mpmath.nstr(worst,4))
# d=5,7 at dps 120
def Aint(p,q):
    return (-1)**p*quad(lambda r: r**p*log(1-mexp(-r))**q/(mexp(r)-1), [0,1,mpmath.inf])
def T_odd(d):
    mid=(d-1)//2
    return sum((-1)**(j+1)*comb(d,j)*Aint(d-j,j) for j in range(mid+1)) - mpf('0.5')*(-1)**(mid+1)*comb(d,mid)*Aint(mid+1,mid)
def piece2(d):
    return mpf('0.5')*quad(lambda w: (-log(w))**d*(1-w)/(1+w)*((1+w)/sqrt(1+w**2)-1)/w, [0,1])
for d in [5,7]:
    mid=(d-1)//2
    T=T_odd(d); p2=piece2(d)
    f1=mpmath.factorial(d)*(1-mpf(2)**(-d))*zeta(d+1)
    tgt=mpmath.factorial(d)/2*E[d+1]
    print(f"d={d}: |T-f1|={mpmath.nstr(abs(T-f1),4)}  |T-p2-tgt|={mpmath.nstr(abs(T-p2-tgt),4)}")
with open('progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: eq(8) dictionary + d=5,7 chain at dps 120\n")

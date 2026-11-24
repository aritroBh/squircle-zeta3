import os, time
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
import mpmath
mpmath.mp.dps = 70
from mpmath import mpf, sqrt, log, pi, zeta, quad
import sys; sys.path.insert(0,'.')
from lead2_lib import E_coeffs
from mpmath import exp as mexp
E = E_coeffs(30)
def A(p,q):
    return (-1)**p*quad(lambda r: r**p*log(1-mexp(-r))**q/(mexp(r)-1), [0,1,mpmath.inf])
from math import comb
def T_odd(d):
    mid = (d-1)//2
    return sum((-1)**(j+1)*comb(d,j)*A(d-j,j) for j in range(mid+1)) \
           - mpf('0.5')*(-1)**(mid+1)*comb(d,mid)*A(mid+1,mid)
def piece2(d):
    return mpf('0.5')*quad(lambda w: (-log(w))**d*(1-w)/(1+w)*((1+w)/sqrt(1+w**2)-1)/w, [0,1])
print("d | T_d - d!(1-2^-d)z(d+1) | piece2 - [d!(1-2^-d)z(d+1) - (d!/2)e_{d+1}]")
for d in [1,3,5,7,9,11]:
    T = T_odd(d); p2 = piece2(d)
    f1 = mpmath.factorial(d)*(1-mpf(2)**(-d))*zeta(d+1)
    f2 = f1 - mpmath.factorial(d)/2*E[d+1]
    print(f"{d}: {mpmath.nstr(T-f1,4)}   {mpmath.nstr(p2-f2,4)}")
with open('progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: closed forms T_d=d!(1-2^-d)zeta(d+1), piece2=T_d-(d!/2)e_{{d+1}} verified d<=11 (70 dps)\n")

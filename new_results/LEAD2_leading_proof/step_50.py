import os, time
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
import mpmath
mpmath.mp.dps = 70
from mpmath import mpf, sqrt, log, pi, zeta, quad
import sys; sys.path.insert(0,'.')
from lead2_lib import E_coeffs, m_e, m_odd, m_oddlast, L1, L2
from math import comb
E = E_coeffs(30)
def Mint(i,j,m): return quad(lambda uu: L1(uu)**i*L2(uu)**j*m(uu), [-1,0,1])
# even d=10
d=10
tot=mpf(0)
for j in range(d//2+1):
    m = m_e if j<d//2 else (lambda uu:(1-uu)/2*m_e(uu))
    tot += (-1)**(j+1)*comb(d,j)*Mint(d-j,j,m)
coll = -mpf('0.5')*quad(lambda uu:(L1(uu)-L2(uu))**d*m_e(uu), [-1,0,1])
tgt = mpmath.factorial(d)/2*E[d+1]
print(f"even d=10: key-coll={mpmath.nstr(tot-coll,4)}  coll-tgt={mpmath.nstr(coll-tgt,4)}")
# odd d=11
d=11; mid=(d-1)//2
tot=mpf(0)
for j in range(mid+1):
    m = m_odd if j<mid else m_oddlast
    tot += (-1)**(j+1)*comb(d,j)*Mint(d-j,j,m)
def Aint(p,q):
    from mpmath import exp as mexp
    return (-1)**p*quad(lambda r: r**p*log(1-mexp(-r))**q/(mexp(r)-1), [0,1,mpmath.inf])
T = sum((-1)**(j+1)*comb(d,j)*Aint(d-j,j) for j in range(mid+1)) - mpf('0.5')*(-1)**(mid+1)*comb(d,mid)*Aint(mid+1,mid)
p2 = mpf('0.5')*quad(lambda w:(-log(w))**d*(1-w)/(1+w)*((1+w)/sqrt(1+w**2)-1)/w, [0,1])
tgt = mpmath.factorial(d)/2*E[d+1]
print(f"odd d=11: T-p2-key={mpmath.nstr(T-p2-tot,4)}  key-tgt={mpmath.nstr(tot-tgt,4)}")
with open('progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: chain checks extended to d=10 (even) and d=11 (odd) at 70 dps\n")

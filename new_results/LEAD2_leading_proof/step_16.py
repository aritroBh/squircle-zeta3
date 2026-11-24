import os, time
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
import mpmath
mpmath.mp.dps = 70
from mpmath import mpf, sqrt, log, pi, zeta, gamma, sin, quad, bernoulli
import sys; sys.path.insert(0, '.')
from lead2_lib import E_coeffs, Efun, m_e, m_odd, m_oddlast, L1, L2
from math import comb

E = E_coeffs(30)
def A(p,q):
    return (-1)**p*quad(lambda r: r**p*log(1-mpmath.exp(-r))**q/(mpmath.exp(r)-1), [0,1,mpmath.inf])
def M(i,j,m):
    return quad(lambda u: L1(u)**i*L2(u)**j*m(u), [-1,0,1])
def T_odd(d):
    mid = (d-1)//2
    return sum((-1)**(j+1)*comb(d,j)*A(d-j,j) for j in range(mid+1)) \
           - mpf('0.5')*(-1)**(mid+1)*comb(d,mid)*A(mid+1,mid)
def piece2(d):
    return mpf('0.5')*quad(lambda w: (-log(w))**d*(1-w)/(1+w)*((1+w)/sqrt(1+w**2)-1)/w, [0,1])
def keysum_odd(d):
    mid = (d-1)//2
    tot = mpf(0)
    for j in range(mid+1):
        m = m_odd if j < mid else m_oddlast
        tot += (-1)**(j+1)*comb(d,j)*M(d-j,j,m)
    return tot
def keysum_even(d):
    tot = mpf(0)
    for j in range(d//2+1):
        m = m_e if j < d//2 else (lambda u: (1-u)/2*m_e(u))
        tot += (-1)**(j+1)*comb(d,j)*M(d-j,j,m)
    return tot
def picsc(d):
    # [z^d] pi/sin(pi z), d odd = 2n-1
    n = (d+1)//2
    return (-1)**(n+1)*(2**(2*n)-2)*bernoulli(2*n)*pi**(2*n)/mpmath.factorial(2*n)
print("=== ODD d: key sum == T_d - piece2 == (d!/2) e_{d+1} ===")
for d in [1,3,5,7,9]:
    ks = keysum_odd(d); T = T_odd(d); p2 = piece2(d)
    tgt = mpmath.factorial(d)/2*E[d+1]
    Tf = mpmath.factorial(d)/2*picsc(d)
    print(f"d={d}: T-p2-key={mpmath.nstr(T-p2-ks,4)}  key-tgt={mpmath.nstr(ks-tgt,4)}  T-Tformula={mpmath.nstr(T-Tf,4)}")
print("=== EVEN d: key sum == -(1/2) int log^d m_e == (d!/2) e_{d+1} ===")
for d in [2,4,6,8]:
    ks = keysum_even(d)
    coll = -mpf('0.5')*quad(lambda u: (L1(u)-L2(u))**d*m_e(u), [-1,0,1])
    tgt = mpmath.factorial(d)/2*E[d+1]
    print(f"d={d}: key-coll={mpmath.nstr(ks-coll,4)}  coll-tgt={mpmath.nstr(coll-tgt,4)}")
with open('progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: odd-d chain d<=9 + even-d collapse d<=8 at 70 dps (T formula via Bernoulli)\n")

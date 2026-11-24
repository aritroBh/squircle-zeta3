import mpmath
mpmath.mp.dps = 80
from mpmath import mpf, sqrt, log, quad, gamma, pi, zeta
from math import comb
lam = log(2)
NN=12
def eta(n): return (1-mpf(2)**(1-n))*zeta(n)
L=[mpf(0)]*(NN+1); L[1]=lam
for n in range(2,NN+1): L[n]=(zeta(n) if n%2==0 else eta(n))/n
E=[mpf(0)]*(NN+1); E[0]=mpf(1)
for n in range(1,NN+1): E[n]=sum(k*L[k]*E[n-k] for k in range(1,n+1))/n
def mainpiece(d):  # -(1/2) ∫_0^1 (-log w)^d (1-w)/sqrt(1+w^2) dw/w
    return -mpf('0.5')*quad(lambda w: (-log(w))**d * (1-w)/sqrt(1+w**2)/w, [0,1])
def boundary(d):
    p = (d+1)//2; q = (d-1)//2
    return (-1)**(p+1) * quad(lambda x: x**p * log(1-mpmath.exp(-x))**q/(mpmath.exp(x)-1), [0,1,mpmath.inf])
for d in [1,3,5,7]:
    tot = mpmath.factorial(d)/2 * E[d+1]
    mp_ = mainpiece(d)
    mid=(d-1)//2
    bterm = -(-1)**(mid+1)*comb(d,mid)/2 * (-1)**((mid+1)+1)*quad(lambda x: x**(mid+1)*log(1-mpmath.exp(-x))**mid/(mpmath.exp(x)-1), [0,1,mpmath.inf])
    # careful: boundary in total = -T_mid(1/(2(1-u))) = -(-1)^{mid+1} binom(d,mid) ∫L1^{mid+1}L2^{mid}/(2(1-u)) du
    # and ∫ L1^{p}L2^{q}/(1-u) du = (-1)^{p+1}∫_0^∞ x^p log^q(1-e^{-x})/(e^x-1) dx, p=mid+1,q=mid
    bterm = -(-1)**(mid+1)*comb(d,mid)/2 * (-1)**(mid+2)*quad(lambda x: x**(mid+1)*log(1-mpmath.exp(-x))**mid/(mpmath.exp(x)-1), [0,1,mpmath.inf])
    print(f"d={d}: total={mpmath.nstr(tot,25)}")
    print(f"   mainpiece={mpmath.nstr(mp_,25)}  total-main={mpmath.nstr(tot-mp_,25)}  bterm={mpmath.nstr(bterm,25)}  sum={mpmath.nstr(mp_+bterm,25)}")

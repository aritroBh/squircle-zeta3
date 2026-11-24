import mpmath
mpmath.mp.dps = 60
from mpmath import mpf, sqrt, log, quad, gamma, pi, zeta
from math import comb
lam = log(2)
NN=14
def eta(n): return (1-mpf(2)**(1-n))*zeta(n)
L=[mpf(0)]*(NN+1); L[1]=lam
for n in range(2,NN+1): L[n]=(zeta(n) if n%2==0 else eta(n))/n
E=[mpf(0)]*(NN+1); E[0]=mpf(1)
for n in range(1,NN+1): E[n]=sum(k*L[k]*E[n-k] for k in range(1,n+1))/n
def piece2(d):  # (1/2)∫_0^inf s^d tanh(s/2)[sqrt2 cosh(s/2)/sqrt(cosh s)-1] ds
    return mpf('0.5')*quad(lambda s: s**d*mpmath.tanh(s/2)*(sqrt(2)*mpmath.cosh(s/2)/sqrt(mpmath.cosh(s))-1), [0,mpmath.inf])
def A(p,q):
    return (-1)**p * quad(lambda r: r**p*log(1-mpmath.exp(-r))**q/(mpmath.exp(r)-1), [0,1,mpmath.inf])
for d in [1,3,5,7]:
    mid = (d-1)//2
    tot = mpmath.factorial(d)/2*E[d+1]
    p1 = sum((-1)**(j+1)*comb(d,j)*A(d-j,j) for j in range(0,mid+1))
    br = mpf('0.5')*(-1)**(mid+1)*comb(d,mid)*A(mid+1,mid)
    p2 = piece2(d)
    print(f"d={d}: total={mpmath.nstr(tot,25)}")
    print(f"   piece1-boundary={mpmath.nstr(p1-br,25)}  piece2={mpmath.nstr(p2,25)}  sum={mpmath.nstr(p1-br+p2,25)}  check={mpmath.nstr(p1-br+p2-tot,8)}")

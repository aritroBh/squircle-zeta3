import mpmath
mpmath.mp.dps = 80
from mpmath import mpf, sqrt, log, quad, gamma, pi, zeta
from math import comb
lam = log(2)
NN=16
def eta(n): return (1-mpf(2)**(1-n))*zeta(n)
L=[mpf(0)]*(NN+1); L[1]=lam
for n in range(2,NN+1): L[n]=(zeta(n) if n%2==0 else eta(n))/n
E=[mpf(0)]*(NN+1); E[0]=mpf(1)
for n in range(1,NN+1): E[n]=sum(k*L[k]*E[n-k] for k in range(1,n+1))/n
# candidate single integrals for odd d
def Jplus(d):   # (1/2)∫(-log w)^d [(1+w)/sqrt(1+w^2)-1] dw/w   (even-d form)
    return mpf('0.5')*quad(lambda w: (-log(w))**d*((1+w)/sqrt(1+w**2)-1)/w, [0,1])
def Jminus(d):  # (1/2)∫(-log w)^d [(1-w)/sqrt(1+w^2)] dw/w -- divergent at 0, skip
    return None
def Jstar(d):   # (1/2)∫(-log w)^d [sqrt(2)cosh(s/2)/sqrt(cosh s)] with w=e^-s: = (1+w)/sqrt(1+w^2): try [(1+w)/sqrt(1+w^2)]/w
    return mpf('0.5')*quad(lambda w: (-log(w))**d*((1+w)/sqrt(1+w**2))/w**1 * 0, [0,1])
for d in [1,3,5,7,9]:
    tot = mpmath.factorial(d)/2*E[d+1]
    jp = Jplus(d)
    print(f"d={d}: total={mpmath.nstr(tot,30)}")
    print(f"   Jplus={mpmath.nstr(jp,30)}  diff={mpmath.nstr(tot-jp,18)}")

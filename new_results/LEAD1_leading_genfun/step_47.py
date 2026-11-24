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
def oddint(d):
    return quad(lambda w: (log(w))**(d+1)*(1+w)/(1+w**2)**mpf('1.5'), [0,1])/(2*(d+1))
for d in [1,3,5,7,9,11]:
    tot = mpmath.factorial(d)/2*E[d+1]
    oi = oddint(d)
    print(f"d={d}: (d!/2)e_{d+1}={mpmath.nstr(tot,30)}  oddint={mpmath.nstr(oi,30)}  diff={mpmath.nstr(tot-oi,8)}")

import mpmath
mpmath.mp.dps = 80
from mpmath import mpf, sqrt, log, quad, gamma, pi, zeta
lam = log(2)
NN = 20
def eta(n): return (1-mpf(2)**(1-n))*zeta(n)
L=[mpf(0)]*(NN+1); L[1]=lam
for n in range(2,NN+1): L[n]=(zeta(n) if n%2==0 else eta(n))/n
E=[mpf(0)]*(NN+1); E[0]=mpf(1)
for n in range(1,NN+1): E[n]=sum(k*L[k]*E[n-k] for k in range(1,n+1))/n
def I(d):
    return quad(lambda w: (-log(w))**d * ((1+w)/sqrt(1+w**2)-1)/w, [0,1])
for d in [0,2,4,6,8,10,12,14,16]:
    val = I(d)
    tgt = mpmath.factorial(d)*E[d+1]
    print(f"d={d}: I_d/(d!)={mpmath.nstr(val/mpmath.factorial(d),30)}  e_{d+1}={mpmath.nstr(E[d+1],30)}  ratio-1={mpmath.nstr(val/tgt-1,6)}")

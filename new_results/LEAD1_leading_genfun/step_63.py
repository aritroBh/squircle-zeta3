import json, mpmath, time
os_nice = __import__('os').nice(10)
mpmath.mp.dps = 250
from mpmath import mpf, sqrt, log, quad, factorial
from fractions import Fraction as Fr
a = 1+sqrt(2); b = sqrt(2)-1
POLES = {'-a': -a, '0': mpf(0), '-b': -b, 'a': a}
D9 = json.load(open('Dtop_d9.json'))
def tval(v): return (v-b)*(v+a)/(4*v)
def uval(v): return (a-v)*(v+b)/(4*v)
def integrand(v, i, j, r):
    t = tval(v); u = uval(v)
    return log(t)**i * log(u)**j / (v - r)
tot = mpf(0); t0=time.time()
mids = [b, mpf('0.5'), 1, mpf('1.7'), a]
for key, rec in D9.items():
    i,j = map(int,key.split(','))
    for pole, pr in rec.items():
        cf = Fr(pr['Dtop'][0])
        if cf:
            val = quad(lambda v: integrand(v,i,j,POLES[pole]), mids)
            val = mpmath.re(val)
            tot += mpf(cf.numerator)/cf.denominator * val
c10 = tot/factorial(8)
lam = log(2)
from mpmath import zeta, pi
NN=14
def eta(n): return (1-mpf(2)**(1-n))*zeta(n)
L=[mpf(0)]*(NN+1); L[1]=lam
for n in range(2,NN+1): L[n]=(zeta(n) if n%2==0 else eta(n))/n
E=[mpf(0)]*(NN+1); E[0]=mpf(1)
for n in range(1,NN+1): E[n]=sum(k*L[k]*E[n-k] for k in range(1,n+1))/n
law = mpf(7936)/(2*factorial(8))*E[10]
print('c(10) =', mpmath.nstr(c10, 50))
print('law   =', mpmath.nstr(law, 50))
print('diff  =', mpmath.nstr(c10-law, 8), f'({time.time()-t0:.0f}s)')

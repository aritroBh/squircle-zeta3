import os, json, time, mpmath
os.nice(10)
mpmath.mp.dps = 300
from mpmath import mpf, sqrt, log, quad, pi, zeta, factorial
from fractions import Fraction as Fr
os.chdir('<project>/research/08_theory/LEAD1_leading_genfun')
a = 1+sqrt(2); b = sqrt(2)-1
POLES = {'-a': -a, '0': mpf(0), '-b': -b, 'a': a}
D9 = json.load(open('Dtop_d9.json'))
def integrand(v, i, j, r):
    t = (v-b)*(v+a)/(4*v); u = (a-v)*(v+b)/(4*v)
    return log(t)**i * log(u)**j / (v - r)
tot = mpf(0); t0=time.time()
for key, rec in D9.items():
    i,j = map(int,key.split(','))
    for pole, pr in rec.items():
        cf = Fr(pr['Dtop'][0])
        if cf:
            val = quad(lambda v: integrand(v,i,j,POLES[pole]), [b,1,a])
            tot += mpf(cf.numerator)/cf.denominator * val
    print(f"key {key} done {time.time()-t0:.0f}s", flush=True)
c10 = tot/factorial(8)
print('c(10) =', mpmath.nstr(c10, 60))
# law
lam = log(2)
NN=14
def eta(n): return (1-mpf(2)**(1-n))*zeta(n)
L=[mpf(0)]*(NN+1); L[1]=lam
for n in range(2,NN+1): L[n]=(zeta(n) if n%2==0 else eta(n))/n
E=[mpf(0)]*(NN+1); E[0]=mpf(1)
for n in range(1,NN+1): E[n]=sum(k*L[k]*E[n-k] for k in range(1,n+1))/n
law = (-1)**(9//2)*7936/(2*factorial(8))*E[10]  # eps_8=+1? eps_8 = [1,-1,-1,1][8%4]=+1; A_9 = 7936
law = Fr(1)*7936/(2*factorial(8))*E[10]
law = mpf(7936)/(2*factorial(8))*E[10]
print('law   =', mpmath.nstr(law, 60))
print('diff  =', mpmath.nstr(c10-law, 8))
json.dump({'c10': mpmath.nstr(c10, 290)}, open('c10_dps300.json','w'))
with open('progress.log','a') as f:
    f.write(f'step_58: c(10) from exact D_9 table vs law: diff printed above ({time.time()-t0:.0f}s)\n')

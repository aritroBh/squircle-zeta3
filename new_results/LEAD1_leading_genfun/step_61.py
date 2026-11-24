import mpmath
mpmath.mp.dps = 300
from mpmath import mpf, log, pi, zeta, factorial
lam = log(2)
NN=14
def eta(n): return (1-mpf(2)**(1-n))*zeta(n)
L=[mpf(0)]*(NN+1); L[1]=lam
for n in range(2,NN+1): L[n]=(zeta(n) if n%2==0 else eta(n))/n
E=[mpf(0)]*(NN+1); E[0]=mpf(1)
for n in range(1,NN+1): E[n]=sum(k*L[k]*E[n-k] for k in range(1,n+1))/n
law = mpf(7936)/(2*factorial(8))*E[10]
c10 = mpf('0.0984598871418917972706502940864227127442040035844511774009152')
print('law   =', mpmath.nstr(law, 60))
print('diff  =', mpmath.nstr(c10-law, 8))
with open('progress.log','a') as f:
    f.write('step_59: c(10) = law to '+mpmath.nstr(c10-law,3)+' (dps 300)\n')

from mpmath import *
mp.dps=25
def E(n):
    n=mpf(n); e=1/n
    # u=e^{-s}, v=e^{-t}; sigma=s+t, s in [0,sigma]
    def inner(sig):
        f=lambda s: (sqrt(1+n**2*exp(-2*sig)*(exp(2*e*s)+exp(2*e*(sig-s))))-1)
        return quad(f,[0,sig/2,sig])*exp(-e*sig)
    c=log(n)
    return quad(inner,[0,c-3,c,c+3,c+10,2*c+40,inf])/n**2
for n in [10,20,100,1000,10000,100000]:
    v=E(n); print(n, v, n*v, n*v-sqrt(2))

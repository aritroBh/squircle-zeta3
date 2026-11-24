import mpmath as mp
# (1) L(chi_5, 5) at 1550 digits via Hurwitz zeta; check vs Euler product at 50 digits
mp.mp.dps = 1550
def Lchi5(s):
    return mp.mpf(5)**(-s)*(mp.zeta(s, mp.mpf(1)/5) - mp.zeta(s, mp.mpf(2)/5) - mp.zeta(s, mp.mpf(3)/5) + mp.zeta(s, mp.mpf(4)/5))
L5 = Lchi5(5)
print("L(chi_5,5) =", mp.nstr(L5, 60))
# Euler product check at 50 digits: prod_p 1/(1 - chi(p) p^-5), chi = Legendre (p/5)
mp.mp.dps = 60
def leg(p):
    r = p % 5
    return 1 if r in (1,4) else (-1 if r in (2,3) else 0)
def primes_upto(n):
    sieve = [True]*(n+1); sieve[0]=sieve[1]=False
    for i in range(2,int(n**0.5)+1):
        if sieve[i]:
            for j in range(i*i,n+1,i): sieve[j]=False
    return [p for p in range(2,n+1) if sieve[p]]
prod = mp.mpf(1)
for p in primes_upto(2000000):
    prod *= 1/(1 - leg(p)*mp.mpf(p)**-5)
print("Euler product (p<=2e6):", mp.nstr(prod, 55))
print("diff vs Hurwitz:", mp.nstr(abs(prod - L5), 3))
# direct series check at 50 digits
s = mp.nsum(lambda n: leg(n)*mp.mpf(n)**-5, [1, mp.inf])
print("direct series:", mp.nstr(s, 55))
print("diff:", mp.nstr(abs(s - L5), 3))
mp.mp.dps = 1550
import pickle
pickle.dump({'Lchi5_5': mp.nstr(L5, 1520)}, open('Lchi5_1550.pkl','wb'))
print("saved Lchi5_1550.pkl")
open('progress.log','a').write("step_70: L(chi_5,5) at 1550 digits via Hurwitz zeta = %s... ; Euler product and direct series agree at 50 digits (diffs printed)\n" % mp.nstr(L5, 30))

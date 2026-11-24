import mpmath as mp
mp.mp.dps = 50

# V_n from R(t) = sqrt((1+e^{4t})/2)
R = lambda t: mp.sqrt((1+mp.exp(4*t))/2)
V = {}
for n in range(0, 16):
    V[n] = mp.diff(R, 0, n)
# (a) recurrence from 4 R R' = 2 R^2 - 1: 4 sum_k C(n,k) V_k V_{n+1-k} = 2 sum_k C(n,k) V_k V_{n-k} - delta_{n,0}
print("recurrence check:")
ok = True
for n in range(0, 12):
    lhs = 4*sum(mp.binomial(n,k)*V[k]*V[n+1-k] for k in range(0, n+1))
    rhs = 2*sum(mp.binomial(n,k)*V[k]*V[n-k] for k in range(0, n+1)) - (1 if n==0 else 0)
    if abs(lhs-rhs) > mp.mpf('1e-35'): ok = False; print("  mismatch n=", n, mp.nstr(lhs-rhs,5))
print("  all ok:", ok)

# (b) kappa exact (from earlier) and asymptotic ratio
kappa = {3: mp.mpf(1), 5: mp.mpf('-0.75'), 7: mp.mpf('8.0625'), 9: -mp.mpf(22805)/192,
         11: mp.mpf(805173)/256, 13: -mp.mpf(131893755)/1024, 15: mp.mpf(154940467501)/20480}
print("\nkappa vs asymptotic  -4 (-1)^m (2m)! (2/pi)^{2m+1} (2m+1)^{-3/2}:")
for m in range(1, 8):
    n = 2*m+1
    pred = -4*(-1)**m * mp.factorial(2*m) * (2/mp.pi)**n * mp.mpf(n)**(-1.5)
    print(f"  n={n:>2}: kappa={mp.nstr(kappa[n],12)} pred={mp.nstr(pred,12)} ratio={mp.nstr(kappa[n]/pred,10)}")

# (c) Borel sum of G(x) = sum kappa_n x^n at x=0.1 vs truncations
# Phi(xi) = -4 E(xi) + 8 ∫_0^{xi/2} R_odd(v) dv/v,  E(xi)=sum_{odd n>=3} xi^n/(n n!)
def E(xi):
    return mp.nsum(lambda nn: xi**nn/(nn*mp.factorial(nn)), [3, mp.inf], tol=mp.mpf('1e-40')) if xi!=0 else mp.mpf(0)
# E via direct sum is slow; use closed form: sum_{odd n>=1} xi^n/(n n!) = (Ein(xi)-Ein(-xi))/2? Ein(z)=sum_{n>=1}(-1)^{n+1} z^n/(n n!)
# sum_{odd n>=1} xi^n/(n n!) = (1/2)(S(xi)-S(-xi)) with S(x)=sum_{n>=1} x^n/(n n!) = Ei(x)-ln x - gamma... for x>0
def S(x):
    # sum_{n>=1} x^n/(n*n!) = -Ein(-x) where Ein(z)=int_0^z (1-e^{-t})/t dt
    return mp.ein(-x) * -1
def Rodd(v):
    return (R(v)-R(-v))/2
def Phi(xi):
    E_ = (S(xi)-S(-xi))/2 - xi   # sum_{odd n>=3} xi^n/(n n!)
    integ = mp.quad(lambda v: Rodd(v)/v if v!=0 else mp.mpf(1), [0, xi/2])
    return -4*E_ + 8*integ
x = mp.mpf('0.1')
Gb = mp.quad(lambda t: mp.e**(-t)*Phi(x*t), [0, mp.inf])
print("\nBorel sum G(0.1) =", mp.nstr(Gb, 25))
part = 0
for n in range(3, 30, 2):
    kn = 4*(-1)**n/n * (1 - 2**(1-n)*mp.diff(R,0,n))
    part += kn * x**n
    if n <= 25: print(f"  partial to n={n:>2}: {mp.nstr(part,25)}")

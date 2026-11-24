# step_07.py -- D11c(1)+(2): generating function of kappa_n: V-recurrence, asymptotics, Borel sum
import mpmath as mp
mp.mp.dps = 50
R = lambda t: mp.sqrt((1+mp.exp(4*t))/2)
V = {n: mp.diff(R, 0, n) for n in range(0, 20)}

# (a) recurrence from R R' = 2 R^2 - 1
ok = all(abs(sum(mp.binomial(n,k)*V[k]*V[n+1-k] for k in range(0,n+1))
             - (2*sum(mp.binomial(n,k)*V[k]*V[n-k] for k in range(0,n+1)) - (1 if n==0 else 0))) < mp.mpf('1e-35')
         for n in range(0, 15))
print("(a) recurrence sum_k C(n,k) V_k V_{n+1-k} = 2 sum_k C(n,k) V_k V_{n-k} - delta_{n0}:", ok)

# (b) asymptotics of kappa_n
kappa = {3: mp.mpf(1), 5: mp.mpf('-0.75'), 7: mp.mpf('8.0625'), 9: -mp.mpf(22805)/192,
         11: mp.mpf(805173)/256, 13: -mp.mpf(131893755)/1024, 15: mp.mpf(154940467501)/20480}
print("(b) kappa_{2m+1} ~ -4 (-1)^m (2m)! (2/pi)^{2m+1} (2m+1)^{-3/2}:")
for m in range(1, 8):
    n = 2*m+1
    pred = -4*(-1)**m*mp.factorial(2*m)*(2/mp.pi)**n*mp.mpf(n)**(-1.5)
    print(f"  n={n:>2}: kappa={mp.nstr(kappa[n],12)} pred={mp.nstr(pred,12)} ratio={mp.nstr(kappa[n]/pred,8)}")

# (c) Borel sum of G(x)=sum kappa_n x^n at x=0.1 vs optimal truncation
def S(x): return mp.quad(lambda u: (mp.exp(u)-1)/u if u!=0 else mp.mpf(1), [0, x])
def Rodd(v): return (R(v)-R(-v))/2
def Phi(xi):
    E_ = (S(xi)-S(-xi))/2 - xi
    integ = mp.quad(lambda v: Rodd(v)/v if v!=0 else mp.mpf(1), [0, xi/2])
    return -4*E_ - 4*xi + 8*integ
x = mp.mpf('0.1')
Gb = mp.quad(lambda t: mp.e**(-t)*Phi(x*t), [0, mp.inf])
part = 0; best = None; bestterm = mp.mpf('1e9')
for n in range(3, 60, 2):
    kn = 4*(-1)**n/n*(1-2**(1-n)*mp.diff(R,0,n)); part += kn*x**n
    if abs(kn*x**n) < bestterm: bestterm = abs(kn*x**n); best = (n, part)
print("\n(c) Borel sum G(0.1) =", mp.nstr(Gb, 25))
print("    optimal truncation n=", best[0], ":", mp.nstr(best[1], 25), " least term", mp.nstr(bestterm, 6))

# (d) closed GF identity at small x: G(x) = 4x - 4 atanh x + 8 sum_{odd n>=3} V_n (x/2)^n / n
x2 = mp.mpf('0.001')
direct = 4*x2 - 4*mp.atanh(x2) + 8*sum(mp.diff(R,0,n)*(x2/2)**n/n for n in range(3, 30, 2))
via_k = sum((4*(-1)**n/n*(1-2**(1-n)*mp.diff(R,0,n)))*x2**n for n in range(3, 30, 2))
print("\n(d) G(0.001) direct:", mp.nstr(direct, 30))
print("    G(0.001) via kappa series:", mp.nstr(via_k, 30))

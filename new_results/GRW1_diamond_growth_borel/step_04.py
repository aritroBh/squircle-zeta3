
import mpmath as mp
mp.mp.dps = 60
from fractions import Fraction
from math import factorial
import json

# ---------- exact rational H_n (square), N=200 ----------
N = 200
KK = N//2+2
a = [Fraction(1, factorial(2*k)) for k in range(KK)]
s = [Fraction(0)]*KK; s[0] = Fraction(1)
for k in range(1, KK):
    s[k] = (a[k] - sum(s[i]*s[k-i] for i in range(1, k)))/2
Hf = [Fraction(0)]*(N+1)
for n in range(N+1):
    Hf[n] = 2*sum(s[k]*Fraction((-1)**(n-2*k), 2**(n-2*k)*factorial(n-2*k)) for k in range(n//2+1))
H = [mp.mpf(v.numerator)/mp.mpf(v.denominator) for v in Hf]

# ---------- f_n from D11b2 450-digit table ----------
d = json.load(open("<project>/research/08_theory/D11b2_stokes_jump/step2_num450_N180.json"))
f450 = [mp.mpf(v) for v in d['f450']]

# ---------- M_{n,j} by quadrature: int_0^inf e^{-u} u^{n-j} k_j(u) du ----------
def kj(j, u):
    if u == 0: return (mp.mpf('0.5'))**2 * (mp.log(2))**j / mp.factorial(j) * ((-1)**j)
    L = -mp.log(1+mp.exp(-u))
    return (1+mp.exp(-u))**(-2) * L**j / mp.factorial(j)

def M(n, j):
    # int_0^inf e^{-u} u^{n-j} k_j(u) du ; peak at u=n-j
    fj = lambda u: mp.exp(-u) * u**(n-j) * kj(j, u)
    pts = [0, max(1,(n-j)//2), max(1,n-j), 2*(n-j)+10, 4*(n-j)+40]
    return mp.quad(fj, pts)

print("=== check 1: M_{n,0} = n! * eta(n) (Dirichlet eta) ===")
for n in [2, 5, 10, 40]:
    lhs = M(n,0); rhs = mp.factorial(n)*(1-mp.mpf(2)**(1-n))*mp.zeta(n)
    print(f"n={n}: M={mp.nstr(lhs,25)}  n!eta(n)={mp.nstr(rhs,25)}  rel={mp.nstr(abs(lhs-rhs)/abs(rhs),3)}")

print()
print("=== check 2: tail bound |M_{n,j}| <= (n-j)!/(j! (j+1)^{n-j+1}) ===")
worst = mp.mpf(0); worstpair=None
for n in [5, 10, 20, 40]:
    for j in range(1, n+1):
        mv = M(n,j)
        bd = mp.factorial(n-j)/(mp.factorial(j)*(j+1)**(n-j+1))
        ratio = abs(mv)/bd
        if ratio > worst: worst, worstpair = ratio, (n,j)
print(f"worst |M|/bound over tested grid = {mp.nstr(worst,6)} at {worstpair}  (must be <=1)")

print()
print("=== check 3: f_n = sum_{j=0}^n H_{n-j} M_{n,j} ===")
for n in [1, 2, 3, 5, 8, 13, 21, 34]:
    tot = mp.mpf(0)
    for j in range(0, n+1):
        tot += H[n-j]*M(n,j)
    fv = f450[n]
    print(f"n={n}: sum={mp.nstr(tot,25)}")
    print(f"      f_n={mp.nstr(fv,25)}   rel diff={mp.nstr(abs(tot-fv)/abs(fv),3)}")

print()
print("=== check 4: relative size of K-tail vs main term ===")
for n in [10, 20, 40, 80, 120, 180]:
    main = H[n]*mp.factorial(n)
    # tail upper bound (proven): sum_{j>=1} |H_{n-j}| (n-j)!/(j!(j+1)^{n-j+1})
    tl = mp.mpf(0)
    for j in range(1, n+1):
        tl += abs(H[n-j])*mp.factorial(n-j)/(mp.factorial(j)*(j+1)**(n-j+1))
    print(f"n={n}: |main|={mp.nstr(abs(main),6)}  tailbound={mp.nstr(tl,6)}  ratio={mp.nstr(tl/abs(main),6)}")

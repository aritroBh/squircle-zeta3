import mpmath as mp, json
mp.mp.dps = 100
from fractions import Fraction
from math import factorial, comb
N=60
KK=N+2
a=[Fraction(0)]*KK; a[0]=Fraction(1)
for k in range(1,KK): a[k]=Fraction((-2)**k,2*factorial(k))
s=[Fraction(0)]*KK; s[0]=Fraction(1)
for k in range(1,KK): s[k]=(a[k]-sum(s[i]*s[k-i] for i in range(1,k)))/2
Hmp=[mp.mpf((2*s[n]).numerator)/mp.mpf((2*s[n]).denominator) for n in range(N+1)]
f450=[mp.mpf(v) for v in json.load(open('../D11b2_stokes_jump/step2_num450_N180.json'))['f450']]

def kjfun(j,u):
    if u==0: return mp.mpf('0.25')*(-mp.log(2))**j/mp.factorial(j)
    L=-mp.log(1+mp.e**(-u))
    return (1+mp.e**(-u))**(-2)*L**j/mp.factorial(j)
def M(n,j):
    f=lambda u: mp.e**(-u)*u**(n-j)*kjfun(j,u)
    pts=[0, max(1,(n-j)//3), max(1,n-j), 2*max(1,n-j)+10]
    return mp.quad(f, pts)
print("=== M_{n,0} = n! eta_D(n) ===")
for n in [2,7,15,40]:
    lhs=M(n,0); rhs=mp.factorial(n)*(1-mp.mpf(2)**(1-n))*mp.zeta(n)
    print(f"n={n}: rel diff={mp.nstr(abs(lhs-rhs)/rhs,3)}")
print("=== |M_{n,j}| <= (n-j)!/(j!(j+1)^{n-j+1}) : worst ratio ===")
worst=mp.mpf(0); wp=None
for n in [6,12,25,50]:
    for j in range(1,n+1):
        mv=M(n,j); bd=mp.factorial(n-j)/(mp.factorial(j)*(j+1)**(n-j+1))
        if abs(mv)/bd>worst: worst, wp = abs(mv)/bd, (n,j)
print("worst ratio:", mp.nstr(worst,6), "at", wp)
print("=== Lemma 2: f_n = sum_j H_{n-j} M_{n,j} ===")
for n in [1,2,3,7,12,20,34,55]:
    tot=mp.mpf(0)
    for j in range(0,n+1):
        tot += Hmp[n-j]*M(n,j)
    print(f"n={n}: sum={mp.nstr(tot,16)} f_n={mp.nstr(f450[n],16)} rel={mp.nstr(abs(tot-f450[n])/abs(f450[n]),3)}")
with open('progress.log','a') as f:
    f.write("step_22: independent Lemma 2 + Lemma 1(b) verification (own quadrature), n up to 55\n")

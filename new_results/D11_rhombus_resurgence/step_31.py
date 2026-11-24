# Fixed composition sign; rerun full independent check
import sympy as sp, json
from mpmath import mp, mpf, quad, log
from math import comb
u, lt, l1 = sp.symbols('u lt l1')
NN = 10
def T(poly, N):
    p = sp.expand(poly); return sum(p.coeff(u, n)*u**n for n in range(0, N+1))
X = T((1+u)*sp.exp(u*lt).series(u, 0, NN+1).removeO(), NN)
Y = T((1+u)*sp.exp(u*l1).series(u, 0, NN+1).removeO(), NN)
W = T(X*X + Y*Y, NN)
Sc = [sp.sqrt(W.coeff(u, 0))]
for n in range(1, NN+1):
    Sc.append(sp.expand((W.coeff(u, n) - sum(Sc[i]*Sc[n-i] for i in range(1, n)))/(2*Sc[0])))
scr = [sp.expand(sc/sp.sqrt(2)) for sc in Sc]
mp.dps = 90
I = []
for n in range(0, NN+1):
    f = sp.lambdify((lt, l1), scr[n], 'mpmath')
    I.append(quad(lambda t: f(log(t), log(1-t)), [0, mpf(1)/32, mpf(1)/4, mpf(1)/2]))
a_dir = [8*I[0]] + [8*sum(I[m]*((-1)**n)*comb(n-1, m-1) for m in range(1, n+1)) for n in range(1, NN+1)]
num = json.load(open('step2_num.json'))
print("INDEPENDENT arclength-integrand-expansion check vs exact closed forms:")
for n in range(0, NN+1):
    an = mpf(num['a'][n]); rd = abs(a_dir[n]-an)/max(1, abs(an))
    print(f"a_{n:2d}  {mp.nstr(a_dir[n], 42)}   rel diff {mp.nstr(rd, 3)}")
# also print I_2, I_3 in closed form terms
print("\nI_2 =", mp.nstr(I[2], 40), "  pi^2/48 =", mp.nstr(mp.pi**2/48, 40))
print("I_3 =", mp.nstr(I[3], 40), "  -zeta(3)/8 =", mp.nstr(-mp.zeta(3)/8, 40))

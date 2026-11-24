import sympy as sp, pickle, time
t0=time.time()
work = '<project>/research/08_theory/TOW5_tower_depth10'
import sys; sys.path.insert(0, work)
d, m = sp.symbols('d m')
i = sp.I
raw = pickle.load(open(work+'/R5_raw.pkl','rb'))
num5 = sp.sympify(raw['num5'], locals={'d':d,'m':m,'I':i})
# Is num5(1,m), num5(2,m)... identically zero in m?
for d0 in [1,2,3,4]:
    v = sp.expand(num5.subs(d, d0))
    print(f"num5({d0},m) == 0 identically in m: {v==0}", f"(terms: {len(sp.Poly(v,m).terms()) if v!=0 else 0})")
# gcd check
den5 = sp.sympify(raw['den5'], locals={'d':d,'m':m,'I':i})
g = sp.gcd(sp.Poly(num5, d, m).as_expr(), sp.Poly(den5, d).as_expr())
print("gcd(num5, den5) =", sp.factor(g))

# Beefed-up direct check at larger d, all keys, plus off-range m
from math import comb, factorial
def stir2(n, k):
    if n == k == 0: return 1
    if n <= 0 or k <= 0 or k > n: return 0
    return sum((-1)**(k-j) * comb(k, j) * j**n for j in range(k+1)) // factorial(k)
def S_poly(dd, pp):
    qq = dd - pp
    t = sp.Symbol('t')
    W = t**2 + (1-t)**2
    out = sp.Integer(0)
    for j in range(1, dd+1):
        Cj = sp.Rational(sp.binomial(sp.Rational(1,2), j) * factorial(j)) * 2**dd
        slo = max(0, j-qq); shi = min(j, pp)
        for s in range(slo, shi+1):
            Sps = stir2(pp, s); Sqj = stir2(qq, j-s)
            if Sps == 0 or Sqj == 0: continue
            out += Cj * Sps * Sqj * t**(2*(j-s)) * (1-t)**(2*s) * W**(dd-j)
    return sp.expand(out) / (factorial(pp)*factorial(qq))
tstar = (1+i)/2
R5 = sp.cancel(num5/den5)
match=0; tot=0
for dd in [12, 15, 18, 22, 26, 30]:
    for pp in range((dd+1)//2, dd+1):
        Sp = S_poly(dd, pp)
        direct = sp.simplify(sp.diff(Sp,'t',5).subs('t',tstar) / Sp.subs('t',tstar))
        pred = sp.simplify(R5.subs({d:dd,m:pp}))
        tot+=1
        if sp.simplify(direct-pred)==0: match+=1
        else: print("MISM", dd, pp)
print(f"beefed direct check d in {{12,15,18,22,26,30}} all keys: {match}/{tot} match ({round(time.time()-t0,1)}s)")

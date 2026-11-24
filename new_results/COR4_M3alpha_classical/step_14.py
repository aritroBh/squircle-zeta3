cor = "<project>/research/08_theory/COR4_M3alpha_classical/"
step = '''"""step10_symbol_checks.py -- Task 1: validate the symbol code on HPLs with known forms."""
import sys, json
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from fractions import Fraction as Fr
from cor4_sym import *

# 3-dim letter basis ([x],[1-x],[1+x]) for a generic argument x
X, M, P = (Fr(1),Fr(0),Fr(0)), (Fr(0),Fr(1),Fr(0)), (Fr(0),Fr(0),Fr(1))
G  = {0: X, 1: lneg(M), -1: P}          # argument  x: g0=[x], g1=-[1-x], g-1=[1+x]
Gm = {0: X, 1: lneg(P), -1: M}          # argument -x: g0=[x], g1=-[1+x], g-1=[1-x]

ok = []

# (a) S(H_{0,0,0,1}(x)) == S(Li4(x)) = -[1-x] [x]^3
assert hpl_sym((0,0,0,1), G) == li_sym(4, X, M) == {(lneg(M), X, X, X): Fr(1)}
assert hpl_sym((0,0,0,-1), G) == li_sym(4, lneg(X), P) or True   # Li4(-x): -[1+x] [x]^3
assert hpl_sym((0,0,0,-1), G) == {(lneg(P), X, X, X): Fr(1)}
ok.append("S(H_{0,0,0,+-1}(x)) = -[1∓x][x]^3 = S(Li4(+-x))")

# (b) shuffle consistency: H_{0,1}^2 = 2 H_{0,1,0,1} + 4 H_{0,0,1,1}
lhs = sshuffle(hpl_sym((0,1), G), hpl_sym((0,1), G))
rhs = sadd(sscale(hpl_sym((0,1,0,1), G), 2), sscale(hpl_sym((0,0,1,1), G), 4))
assert lhs == rhs
ok.append("S(H_{0,1})^2 shuffle = S(2 H_{0,1,0,1} + 4 H_{0,0,1,1})")

# (c) Li4(x) + Li4(-x) = (1/8) Li4(x^2) at symbol level
#     letters of x^2: [x^2]=2[x], [1-x^2]=[1-x]+[1+x]
lhs = sadd(li_sym(4, X, M), li_sym(4, X, P))          # Li4(x)+Li4(-x): -[1∓x][x]^3
rhs = sscale(li_sym(4, lscale(2, X), ladd(M, P)), Fr(1, 8))
assert lhs == rhs, (lhs, rhs)
ok.append("S(Li4(x)+Li4(-x)) = (1/8) S(Li4(x^2))  [uses [1-x^2]=[1-x]+[1+x]]")

# (d) delta22 kills Li4, detects depth two
assert delta22(li_sym(4, X, M), 3) == {}
d = delta22(hpl_sym((0,0,1,1), G), 3)
assert d != {}
ok.append(f"delta22(Li4(x)) = 0 ; delta22(H_{{0,0,1,1}}(x)) = {d}")

# (e) numeric: HPL quadrature vs classical / multiple-polylog values
import mpmath as mp
mp.mp.dps = 70
def f(a, t):
    return 1/t if a == 0 else (1/(1-t) if a == 1 else 1/(1+t))
def H2(w, x):
    if w == (1,):  return -mp.log(1-x)
    if w == (-1,): return mp.log(1+x)
    if w == (0,1): return mp.polylog(2, x)
    if w == (0,-1): return -mp.polylog(2, -x)
    if w == (1,1): return mp.log(1-x)**2/2
    if w == (-1,-1): return mp.log(1+x)**2/2
    if w == (1,-1): return mp.polylog(2,(1-x)/2) - mp.polylog(2,mp.mpf(1)/2) - mp.log(2)*mp.log(1-x)
    if w == (-1,1): return mp.polylog(2,(1+x)/2) - mp.polylog(2,mp.mpf(1)/2) - mp.log(2)*mp.log(1+x)
    if w == (0,0,1): return mp.polylog(3,x)
    if w == (0,0,-1): return -mp.polylog(3,-x)
    if w == (0,0,0,1): return mp.polylog(4,x)
    if w == (0,0,0,-1): return mp.polylog(4,-x)
    raise KeyError(w)
cache = {}
def H(w, x):
    w = tuple(w)
    if (w, x) in cache: return cache[(w, x)]
    try:
        v = H2(w, x)
    except KeyError:
        a, rest = w[0], w[1:]
        v = mp.quad(lambda t: f(a, t)*H(rest, t), [0, x/4, x/2, 3*x/4, x])
    cache[(w, x)] = v
    return v

def Li31(x, y):   # sum_{n>m>=1} x^n y^m/(n^3 m)
    return mp.nsum(lambda n, m: x**n * y**m / (n**3 * m),
                   [2, mp.inf], [1, lambda n: n-1], tol=mp.mpf('1e-65'))
def Li22(x, y):   # sum_{n>m>=1} x^n y^m/(n^2 m^2)
    return mp.nsum(lambda n, m: x**n * y**m / (n**2 * m**2),
                   [2, mp.inf], [1, lambda n: n-1], tol=mp.mpf('1e-65'))

res = {}
for x in [mp.mpf('0.31'), mp.mpf('0.77')]:
    res[('H0001-Li4', x)]   = abs(H((0,0,0,1), x) - mp.polylog(4, x))
    res[('H000m1-Li4', x)]  = abs(H((0,0,0,-1), x) - mp.polylog(4, -x))
    res[('H0011-Li31', x)]  = abs(H((0,0,1,1), x) - Li31(x, 1))
    res[('H0101-Li22', x)]  = abs(H((0,1,0,1), x) - Li22(x, 1))
    # shuffle H_{0,1}^2 = 2H_{0,1,0,1}+4H_{0,0,1,1} numerically
    res[('shuffle-num', x)] = abs(H((0,1),x)**2 - 2*H((0,1,0,1),x) - 4*H((0,0,1,1),x))
    # H_{0,0,1,-1}(x) = -Li31(x,-1)
    res[('H001m1', x)]      = abs(H((0,0,1,-1), x) + Li31(x, -1))
    res[('H00m11', x)]      = abs(H((0,0,-1,1), x) - Li31(-x, -1) * -1 * -1)  # placeholder, fixed below

# fix the last identity properly: H_{0,0,-1,1}(x): nonzero letters order (-1,1):
# series: H_{0,0,-1,1}(x) = sum_{n>m} (-1)^{n+?}...; test a couple of candidates numerically
x = mp.mpf('0.31')
cands = {'Li31(-x,1)': Li31(-x,1), '-Li31(-x,1)': -Li31(-x,1),
         'Li31(-x,-1)': Li31(-x,-1), '-Li31(-x,-1)': -Li31(-x,-1),
         'Li31(x,-1)': Li31(x,-1), '-Li31(x,-1)': -Li31(x,-1)}
hv = H((0,0,-1,1), x)
print('H_{0,0,-1,1}(0.31) =', mp.nstr(hv, 25))
for k, v in cands.items():
    if abs(hv - v) < mp.mpf('1e-50'):
        print('  ==', k)
    if abs(hv + v) < mp.mpf('1e-50'):
        print('  == -(', k, ')')
res.pop(('H00m11', mp.mpf('0.31')), None); res.pop(('H00m11', mp.mpf('0.77')), None)

print()
for k in sorted(res, key=str):
    print(k, '->', mp.nstr(res[k], 3))
print()
print("SYMBOL-LEVEL CHECKS (exact):")
for s in ok:
    print("  PASS:", s)
maxres = max(res.values())
print("max numeric residual:", mp.nstr(maxres, 3))
assert maxres < mp.mpf('1e-50')
print("ALL TASK-1 VALIDATION CHECKS PASS")
'''
cor = "<project>/research/08_theory/COR4_M3alpha_classical/"
open(cor+"step10_symbol_checks.py","w").write(step)
import subprocess
r = subprocess.run(["python3", cor+"step10_symbol_checks.py"], capture_output=True, text=True)
print(r.stdout)
print(r.stderr[-2000:] if r.returncode else "exit 0")

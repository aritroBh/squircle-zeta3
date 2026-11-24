cor = "<project>/research/08_theory/COR4_M3alpha_classical/"
step = '''"""step10_symbol_checks.py -- Task 1: validate the symbol code on HPLs with known forms."""
import sys, json
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from fractions import Fraction as Fr
from cor4_sym import *
from hpl_series import H_series, Li31, Li22

# 3-dim letter basis ([x],[1-x],[1+x]) for a generic argument x
X, M, P = (Fr(1),Fr(0),Fr(0)), (Fr(0),Fr(1),Fr(0)), (Fr(0),Fr(0),Fr(1))
G  = {0: X, 1: lneg(M), -1: P}          # argument  x
Gm = {0: X, 1: lneg(P), -1: M}          # argument -x

ok = []

# (a) S(H_{0,0,0,1}(x)) == S(Li4(x)) = -[1-x] [x]^3 ; H_{0,0,0,-1}(x) = -Li4(-x)
assert canon(hpl_sym((0,0,0,1), G), 3) == canon(li_sym(4, X, M), 3) == {(1, 0, 0, 0): Fr(-1)}
assert canon(hpl_sym((0,0,0,-1), G), 3) == {(2, 0, 0, 0): Fr(1)}
assert canon(hpl_sym((0,0,0,-1), G), 3) == canon(sscale(li_sym(4, X, P), -1), 3)
ok.append("S(H_{0,0,0,1}) = -[1-x][x]^3 = S(Li4(x));  S(H_{0,0,0,-1}) = +[1+x][x]^3 = -S(Li4(-x))")

# (b) shuffle consistency: H_{0,1}^2 = 2 H_{0,1,0,1} + 4 H_{0,0,1,1}
lhs = canon(sshuffle(hpl_sym((0,1), G), hpl_sym((0,1), G)), 3)
rhs = canon(sadd(sscale(hpl_sym((0,1,0,1), G), 2), sscale(hpl_sym((0,0,1,1), G), 4)), 3)
assert lhs == rhs
ok.append("S(H_{0,1})^2 shuffle = S(2 H_{0,1,0,1} + 4 H_{0,0,1,1})")

# (c) Li4(x) + Li4(-x) = (1/8) Li4(x^2) at symbol level
lhs = canon(sadd(li_sym(4, X, M), li_sym(4, X, P)), 3)
rhs = canon(sscale(li_sym(4, lscale(2, X), ladd(M, P)), Fr(1, 8)), 3)
assert lhs == rhs
ok.append("S(Li4(x)+Li4(-x)) = (1/8) S(Li4(x^2))  [uses [1-x^2]=[1-x]+[1+x]]")

# (d) delta22: kills Li4 and H_{0,0,1,1}=Li31(x,1) (which is only seen by the (1,3) part),
#     detects H_{0,1,0,1}=Li22(x,1)
assert delta22(li_sym(4, X, M), 3) == {}
assert delta22(hpl_sym((0,0,1,1), G), 3) == {}
d = delta22(hpl_sym((0,1,0,1), G), 3)
assert d == {((0, 1), (0, 1)): Fr(1)}, d
ok.append("delta22(Li4(x)) = 0; delta22(Li31(x,1)) = 0; delta22(Li22(x,1)) = ((01),(01)) != 0")

# (e) numeric: series HPL vs classical / multiple-polylog values (dps 60)
import mpmath as mp
mp.mp.dps = 60
res = {}
for x in [mp.mpf('0.31'), mp.mpf('0.77')]:
    res['H0001-Li4']  = max(res.get('H0001-Li4',0), abs(H_series((0,0,0,1), x) - mp.polylog(4, x)))
    res['H000m1+Li4'] = max(res.get('H000m1+Li4',0), abs(H_series((0,0,0,-1), x) + mp.polylog(4, -x)))
    res['H0011=Li31(x,1)']  = max(res.get('H0011=Li31(x,1)',0), abs(H_series((0,0,1,1), x) - Li31(x, 1)))
    res['H0101=Li22(x,1)']  = max(res.get('H0101=Li22(x,1)',0), abs(H_series((0,1,0,1), x) - Li22(x, 1)))
    res['H001m1=-Li31(x,-1)'] = max(res.get('H001m1=-Li31(x,-1)',0), abs(H_series((0,0,1,-1), x) + Li31(x, -1)))
    res['shuffle-num'] = max(res.get('shuffle-num',0),
        abs(H_series((0,1),x)**2 - 2*H_series((0,1,0,1),x) - 4*H_series((0,0,1,1),x)))
# identify H_{0,0,-1,1}(x) and H_{0,-1,0,1}(x)
x = mp.mpf('0.31')
hv = H_series((0,0,-1,1), x)
for name, v in [('Li31(-x,-1)', Li31(-x,-1)), ('-Li31(-x,-1)', -Li31(-x,-1)),
                ('Li31(-x,1)', Li31(-x,1)), ('-Li31(-x,1)', -Li31(-x,1))]:
    if abs(hv - v) < mp.mpf('1e-50'):
        print(f'H_{{0,0,-1,1}}(x) == {name}   (residual {mp.nstr(abs(hv-v),3)})')
hv2 = H_series((0,-1,0,1), x)
for name, v in [('Li22(-x,-1)', Li22(-x,-1)), ('-Li22(-x,-1)', -Li22(-x,-1)),
                ('Li22(-x,1)', Li22(-x,1)), ('-Li22(-x,1)', -Li22(-x,1))]:
    if abs(hv2 - v) < mp.mpf('1e-50'):
        print(f'H_{{0,-1,0,1}}(x) == {name}   (residual {mp.nstr(abs(hv2-v),3)})')

print()
print("numeric residuals (should be < 1e-55):")
for k, v in sorted(res.items()):
    print(f"  {k:>20}: {mp.nstr(v, 3)}")
print()
print("SYMBOL-LEVEL CHECKS (exact):")
for s in ok:
    print("  PASS:", s)
assert max(res.values()) < mp.mpf('1e-55')
print("ALL TASK-1 VALIDATION CHECKS PASS")
'''
open(cor+"step10_symbol_checks.py","w").write(step)

# hpl_series.py module (shared)
hplmod = '''"""hpl_series.py -- fast high-precision HPL evaluation via exact recursive series,
and multiple polylogs Li_{3,1}(x,y)=sum_{n>m} x^n y^m/(n^3 m), Li_{2,2} likewise.
Cross-checked against the quadrature of G1_zeta3_any_angle/hpl_eval.py (step15)."""
import math
import mpmath as mp

def hpl_coeffs(w, N):
    """coefficients c[0..N] of H_w(x) = sum_{n>=1} c_n x^n (c_0 = 0 for nonempty w)."""
    if len(w) == 0:
        c = [mp.mpf(0)] * (N + 1); c[0] = mp.mpf(1); return c
    a, rest = w[0], w[1:]
    c = hpl_coeffs(rest, N)
    out = [mp.mpf(0)] * (N + 1)
    if a == 0:
        for n in range(1, N + 1):
            out[n] = c[n] / n
    elif a == 1:
        S = mp.mpf(0)
        for n in range(1, N + 1):
            S += c[n - 1]
            out[n] = S / n
    else:
        A = mp.mpf(0)
        for n in range(1, N + 1):
            A = c[n - 1] - A
            out[n] = A / n
    return out

def H_series(w, x, N=None):
    w = tuple(w)
    if w == (0,):
        return mp.log(x)
    if N is None:
        N = int((mp.mp.dps + 15) / abs(math.log10(abs(x)))) + 10
    c = hpl_coeffs(w, N)
    tot = mp.mpf(0); xn = x
    for n in range(1, N + 1):
        tot += c[n] * xn; xn *= x
    return tot

def _mpl_series(x, y, a, b):
    tot = mp.mpf(0); inner = mp.mpf(0); xn = x; n = 1
    while True:
        if n >= 2:
            m = n - 1
            inner += y**m / mp.mpf(m)**b
        term = xn * inner / mp.mpf(n)**a
        tot += term
        if n > 40 and abs(term) < mp.mpf(10) ** (-mp.mp.dps + 10):
            break
        n += 1; xn *= x
    return tot

def Li31(x, y): return _mpl_series(x, y, 3, 1)
def Li22(x, y): return _mpl_series(x, y, 2, 2)
'''
open(cor+"hpl_series.py","w").write(hplmod)
import sys, subprocess
r = subprocess.run([sys.executable, cor+"step10_symbol_checks.py"], capture_output=True, text=True, timeout=600)
print(r.stdout)
print(r.stderr[-1500:] if r.returncode else "exit 0")

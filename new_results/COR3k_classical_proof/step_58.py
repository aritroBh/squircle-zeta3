import mpmath as mp, ast, time
from cor3_common import *
import sympy as sp

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
CD = load_classical_dict(G1 + "classical_dict.json")
mp.mp.dps = 65
lam = mp.log(2)

def f_at(a, t):
    return 1/t if a == 0 else (1/(1-t) if a == 1 else 1/(1+t))
def C2num(ab, t):   # proved weight-2 classical forms (mpmath)
    if ab == (0, 1): return mp.polylog(2, t)
    if ab == (0, -1): return -mp.polylog(2, -t)
    if ab == (1, -1): return -mp.polylog(2, (1+t)/2) - mp.log(1-t)*mp.log(1+t) + lam*mp.log(1+t) - lam**2/2 + mp.pi**2/12
    if ab == (-1, 1): return mp.polylog(2, (1+t)/2) - lam*mp.log(1+t) + lam**2/2 - mp.pi**2/12
def Hw3_quad(w, X):  # defining integral, inner weight-2 via proved form
    a1, a2, a3 = w
    return mp.quad(lambda t: f_at(a1, t)*C2num((a2, a3), t), [0, X])
def Hw3_nested(w, X):  # fully nested defining integral (cross-check)
    a1, a2, a3 = w
    H1 = lambda s: (mp.log(s) if a3 == 0 else (-mp.log(1-s) if a3 == 1 else mp.log(1+s)))
    return mp.quad(lambda t: f_at(a1, t)*mp.quad(lambda s: f_at(a2, s)*H1(s), [0, t]), [0, X])

# cross-check the two quadratures once
d = abs(Hw3_quad((0,1,-1), mp.mpf('0.7')) - Hw3_nested((0,1,-1), mp.mpf('0.7')))
print("quad cross-check (proved weight-2 vs fully nested):", mp.nstr(d, 3))

print("\n=== NUMERICAL (60-digit): H_w quad vs C_w ===")
for w_str in sorted(CD.keys()):
    w = ast.literal_eval(w_str)
    Cf = sp.lambdify(x, CD[w_str], 'mpmath')
    for X in [mp.mpf('0.3'), mp.mpf('0.7'), mp.mpf('0.95')]:
        hv = Hw3_quad(w, X); cv = Cf(X)
        print(f"  {w_str:>12} x={X}: |quad - C_w| = {mp.nstr(abs(hv-cv), 3)}")

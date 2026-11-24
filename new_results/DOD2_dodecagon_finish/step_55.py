code = r'''
import mpmath as mp, pickle, time
import sympy as sp
from sympy import Rational, sqrt, symbols, expand, Function, simplify
mp.mp.dps = 620
t0=time.time()
D1 = '../DOD1_dodecagon/'
r3 = sqrt(3)
l2,l3,lu = symbols('l2 l3 lu')
Li2A = Function('Li2A'); Li3A = Function('Li3A')
# exact x-points
Xs = {1: 1-2*r3/3, 2: 1-r3, 5: 1-r3/3, 6: 1-r3/2}
# ---- exact log factorization: theta in Q(sqrt3), theta>0 -> (a,b,m) with theta = 2^a 3^b (2+sqrt3)^m
UP = 2+r3  # fundamental unit >1
def logfac(theta):
    th = sp.simplify(sp.expand(theta))
    assert th.is_real
    sgn = 1
    if th < 0: sgn = -1; th = -th
    # numerical solve for (a,b,m)
    lt = sp.log(th)
    import mpmath as mpx
    mpx.mp.dps = 80
    lv = mpx.mpf(sp.N(lt, 70)); L2 = mpx.log(2); L3 = mpx.log(3); LU = mpx.log(2-mpx.sqrt(3))
    rel = mpx.pslq([lv, L2, L3, LU], tol=mpx.mpf('1e-70'), maxcoeff=10**6, maxsteps=2000)
    assert rel is not None and rel[0] != 0
    a = -Rational(rel[1], rel[0]); b = -Rational(rel[2], rel[0]); c = -Rational(rel[3], rel[0])
    # c is coefficient of lu = log(2-sqrt3) = -log(2+sqrt3); verify exactly
    cand = 2**a * 3**b * (2-r3)**c
    ok = sp.simplify(sp.expand(th - cand)) == 0
    if not ok:
        # try with (2+sqrt3) powers via c on lu: lu=log(2-r3) so factor (2-r3)^c
        raise ValueError("logfac failed for %s (cand %s)" % (th, cand))
    return a*l2 + b*l3 + c*lu
# test
for th in [1-Xs[2], abs(1-Xs[5]/Xs[2]), r3, (r3-1), 3-r3, (r3+3)/3]:
    print(th, "->", logfac(th))
print("logfac OK, t=%ds"%(time.time()-t0))
pickle.dump({'done':1}, open('logfac_test.pkl','wb'))
'''
open('step_15.py','w').write(code)
exec(code)

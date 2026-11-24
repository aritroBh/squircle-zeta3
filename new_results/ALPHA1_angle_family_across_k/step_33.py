step04 = r'''"""step_04.py -- component tests: pf_from_poly, pf_fold numerically."""
import sys
sys.path.insert(0, '.')
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 40
from cwfield import make_quad
from cwintegrand import CornerIntegrand
from cwred import Reducer
from cwbuild import build_integrand_pf

c = Fr(1,2); r = Fr(1,3)
F = make_quad(-r)
red = Reducer(F, F(r))
ci = CornerIntegrand(F, F(c), F(r))
rmp = mp.mpf(1)/3
th = mp.mpc(0, 1)/mp.sqrt(3)
def q2mp(q):
    return mp.mpf(int(q.a.numerator))/int(q.a.denominator) + (mp.mpf(int(q.b.numerator))/int(q.b.denominator))*th
def kfun(P):
    def f(w):
        tot = mp.mpc(0)
        for rho, cf in P[0].items():
            rr = mp.mpf(int(rho.a.numerator))/int(rho.a.denominator)
            tot += q2mp(cf)/(w - rr)
        for j, e in enumerate(P[1], start=1):
            tot += q2mp(e)/(w - th)**j
        for j, e in enumerate(P[2], start=1):
            tot += q2mp(e)/(w + th)**j
        tot += q2mp(P[3])
        return tot
    return f

# test 1: pf_from_poly on a k=2 sqrt piece
T = ci.T_poly(2, 2)
num, den = ci.sqrt_piece_nd(2, T)
P = red.pf_from_poly(num, den)
def ndfun(w):
    w = mp.mpc(w)
    nv = sum(q2mp(cf)*w**i for i, cf in enumerate(num))
    dv = sum(q2mp(cf)*w**i for i, cf in enumerate(den))
    return nv/dv
fP = kfun(P)
for w0 in (mp.mpf('0.4'), mp.mpf('0.7'), mp.mpf('0.9')):
    print(f"pf_from_poly check w={w0}: nd={mp.nstr(ndfun(w0),20)} pf={mp.nstr(fP(w0),20)} diff={mp.nstr(abs(ndfun(w0)-fP(w0)),3)}")

# test 2: pf_fold: F(R)(w) =? r R(r/w)/w^2
fF = kfun(red.pf_fold(P))
for w0 in (mp.mpf('0.4'), mp.mpf('0.7'), mp.mpf('0.9')):
    lhs = fF(w0)
    rhs = rmp * fP(rmp/w0) / w0**2
    print(f"pf_fold check w={w0}: fold={mp.nstr(lhs,20)} rR(r/w)/w^2={mp.nstr(rhs,20)} diff={mp.nstr(abs(lhs-rhs),3)}")
print("STEP04 DONE")
'''
open('step_04.py','w').write(step04)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'step_04.py'], capture_output=True, text=True)
print(r.stdout)
print(r.stderr[-2000:] if r.returncode else "RC=0")

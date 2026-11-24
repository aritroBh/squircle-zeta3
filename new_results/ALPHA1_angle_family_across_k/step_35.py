step06 = r'''"""step_06.py -- test pf_mul_DL and feed-fold numerically at k=2."""
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
om1 = lambda w: 1/(w-rmp) + 1/(w+1) - 1/w
om2 = lambda w: 1/(w-1) + 1/(w+rmp) - 1/w

ipf = build_integrand_pf(ci, red, 2)
for key, W in ipf.items():
    i, j = key
    simple, higher, const = red.pf_split(W)
    G = red.pf_integrate_higher(higher)
    C = red.pf_eval_at(G, red.r)
    Gm = [G[0], G[1], G[2], -C]
    fG = kfun(Gm)
    for DL, om, nm in ((red.DL1, om1, 'DL1'), (red.DL2, om2, 'DL2')):
        feed = red.pf_mul_DL(Gm, DL)
        ff = kfun(feed)
        for w0 in (mp.mpf('0.4'), mp.mpf('0.75')):
            lhs = ff(w0); rhs = fG(w0)*om(w0)
            print(f"key{key} {nm} w={w0}: mulDL={mp.nstr(lhs,15)} G*om={mp.nstr(rhs,15)} diff={mp.nstr(abs(lhs-rhs),2)}")
        # fold of feed
        ffold = kfun(red.pf_fold(feed))
        for w0 in (mp.mpf('0.4'), mp.mpf('0.75')):
            lhs = ffold(w0); rhs = rmp*ff(rmp/w0)/w0**2
            print(f"key{key} {nm} fold w={w0}: fold={mp.nstr(lhs,15)} rF(r/w)/w^2={mp.nstr(rhs,15)} diff={mp.nstr(abs(lhs-rhs),2)}")
print("STEP06 DONE")
'''
open('step_06.py','w').write(step06)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'step_06.py'], capture_output=True, text=True)
print(r.stdout)
print(r.stderr[-2000:] if r.returncode else "RC=0")

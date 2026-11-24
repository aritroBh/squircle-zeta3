import os, subprocess, time
Z3 = 'Zeta' + '(3)'
ROOT = os.path.join('<home>/Downloads/Research', Z3)
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
PY = os.path.join(ROOT, '.venv/bin/python')
step08 = r'''"""step08: INDEPENDENT NUMERIC SANITY (not part of the proof).
G_X = W4pp_X - sum c_k B_k - sum z_j F_j is EXACTLY in the T0 rowspace (stage 2 + verifier).
Evaluate G_rat, G_sq directly with the A3f4 evaluator (N=2000, prec=700, precision set after imports):
both must vanish to ~1e-200, confirming the whole identity independently of the linear algebra."""
import sys, json, pickle, time
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/M4b_w4pp_exact')
import m4b_common as mc
from fractions import Fraction as F
sys.path.insert(0, '<project>/research/08_theory/A3f4_glanois_basis/Q')
import check_conversion as cc
from flint import ctx, acb, arb
import mpmath as mp
ctx.prec = 700; ctx.cap = 2200
cc.roots = [acb(0, arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0], cc.roots[2], cc.roots[4], cc.roots[6] = acb(1), acb(0, 1), acb(-1), acb(0, -1)
cc.N = 2000
mp.mp.dps = 100
def clear():
    for name in ('series', 'half', 'g_at_one', 'h'): getattr(cc, name).cache_clear()
def aval(c):
    c = F(c); return acb(str(c.numerator)) / acb(str(c.denominator))
def vecval(d):
    s = acb(0)
    for w, c in d.items(): s += aval(c) * cc.h(w)
    clear()
    return s
def mag(z):
    return mp.sqrt(mp.mpf(z.real.mid().str(30, radius=False, more=True))**2
                   + mp.mpf(z.imag.mid().str(30, radius=False, more=True))**2)
t0 = time.time()
G = pickle.load(open(str(mc.HERE/'stage2_G.pkl'), 'rb'))
for name in ('Grat', 'Gsq'):
    r = vecval(G[name])
    mc.log(f'{name}: support={len(G[name])} |value|={mp.nstr(mag(r),3)}')
mc.log(f'step08 done s={time.time()-t0:.1f}')
'''
open('step08_numeric_sanity.py','w').write(step08)
t = time.time()
r = subprocess.run([PY, 'step08_numeric_sanity.py'], capture_output=True, text=True, cwd=WORK)
print('wall', round(time.time()-t,1),'s  EXIT', r.returncode)
print(r.stdout)
print(r.stderr[-2000:] if r.returncode else '')

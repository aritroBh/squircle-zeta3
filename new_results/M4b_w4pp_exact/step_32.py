import os
ROOT = '<project>')
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
step03 = '''"""step_03: numeric checks at >>40 digits (sanity only; proof is exact linear algebra).
Trap-safe: precision is set AFTER all imports (check_conversion import caps prec at 200 bits).
Evaluator = A3f4 check_conversion.h, roots recomputed after raising prec; N=2000.
Checks:
  (a) every atom vector B_k vs direct mpmath (log/polylog/zeta products) -- must agree < 1e-45;
  (b) eval(W4pp_rat) + sqrt2*eval(W4pp_sq) vs the M4a G-value (gpl_w4pp_check.json, == quad at 90 dps).
"""
import sys, json, pickle, time
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/M4b_w4pp_exact')
import m4b_common as mc          # binds A3c cmzv/relations into sys.modules first
from fractions import Fraction as F
sys.path.insert(0, '<project>/research/08_theory/A3f4_glanois_basis/Q')
import check_conversion as cc
from flint import ctx, acb, arb
import mpmath as mp
# ---- precision set AFTER all imports (trap) ----
ctx.prec = 700; ctx.cap = 2200
cc.roots = [acb(0, arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0], cc.roots[2], cc.roots[4], cc.roots[6] = acb(1), acb(0, 1), acb(-1), acb(0, -1)
cc.N = 2000
mp.mp.dps = 120
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
W = pickle.load(open(str(mc.HERE/'stage1_w4pp_vectors.pkl'), 'rb'))
AT = pickle.load(open(str(mc.HERE/'stage2_atoms.pkl'), 'rb'))['atoms']

# two-N convergence spot check on 6 words of the target
sw = list(W['rat'])[:3] + list(W['sq'])[:3]
clear(); hi = [cc.h(w) for w in sw]; clear(); cc.N = 1000
lo = [cc.h(w) for w in sw]; clear(); cc.N = 2000
worst = max(mag(a-b) for a, b in zip(hi, lo))
mc.log(f'two-N convergence (6 words, N=2000 vs 1000): max diff={mp.nstr(worst,3)}')

# (a) atoms
from qsqrt2 import QS2 as Q, SQ2, A, B
mq = lambda z: mp.mpf(str(z.p)) + mp.sqrt(2)*mp.mpf(str(z.q))
def atom_value(a):
    if a[0] == 'log': return mp.log(mq(a[1]))
    if a[0] == 'pi2': return mp.pi**2
    if a[0] == 'z3': return mp.zeta(3)
    if a[0] == 'li': return mp.polylog(a[1], mq(a[2]))
    raise ValueError(a)
def decode(enc):
    out = []
    for a in enc:
        if a[0] == 'log': out.append(('log', Q(F(a[1]))))
        elif a[0] in ('pi2', 'z3'): out.append((a[0],))
        elif a[0] == 'li':
            s = a[2]
            import re
            m = re.fullmatch(r"\\((.+?)([+-])(\\d+(?:/\\d+)?)\\*sqrt2\\)", s)
            if m:
                p = F(m.group(1)); q = F(m.group(3)); q = q if m.group(2) == '+' else -q
                out.append(('li', a[1], Q(p, q)))
            elif s.endswith('*sqrt2'):
                out.append(('li', a[1], Q(0, F(s[:-6]))))
            else:
                out.append(('li', a[1], Q(F(s))))
    return out
worst = mp.mpf(0)
for name, enc, v in AT:
    direct = mp.fprod([atom_value(a) for a in decode(enc)])
    got = vecval(v)
    err = abs(mp.mpc(got.real.mid().str(60, radius=False, more=True),
                     got.imag.mid().str(60, radius=False, more=True)) - direct)
    worst = max(worst, err)
    mc.log(f'atom {name}: support={len(v)} err_vs_direct={mp.nstr(err,3)}')
    assert err < mp.mpf("1e-45"), (name, err)
mc.log(f'(a) all {len(AT)} atoms agree with direct mpmath < 1e-45; worst={mp.nstr(worst,3)}')

# (b) W4'' combination
vr = vecval(W['rat']); vs = vecval(W['sq'])
val = mp.mpc(vr.real.mid().str(115, radius=False, more=True), vr.imag.mid().str(115, radius=False, more=True)) \\
    + mp.sqrt(2)*mp.mpc(vs.real.mid().str(115, radius=False, more=True), vs.imag.mid().str(115, radius=False, more=True))
G = mp.mpf(json.load(open(str(mc.THEORY/'M4a_decomposition/gpl_w4pp_check.json')))['G_value'])
err = abs(val - G)
mc.log(f'(b) eval(W4pp_rat)+sqrt2*eval(W4pp_sq) = {mp.nstr(mp.re(val),70)}')
mc.log(f'    imag part = {mp.nstr(mp.im(val),3)}; err vs M4a G-value = {mp.nstr(err,3)}')
assert mp.im(val) < mp.mpf("1e-60") and err < mp.mpf("1e-60")
mc.log(f'step_03 done s={time.time()-t0:.1f}')
'''
with open('step_03.py','w') as f: f.write(step03)
import subprocess, time
PY = '<project>/.venv/bin/python'
t=time.time()
r = subprocess.run([PY, 'step_03.py'], capture_output=True, text=True, cwd=WORK)
print('wall', round(time.time()-t,1),'s')
print(r.stdout)
print(r.stderr[-3000:] if r.returncode else '', 'EXIT', r.returncode)

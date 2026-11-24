base = '<project>'
here = base + '/research/08_theory/M4b2_w4pp_nice'
step09 = '''"""M4b2 step09: 1060-digit numeric sanity (not part of the proof).
Evaluate R''_nice = sum_k (c_rat,k + sqrt2*c_sq,k) B_k directly with mpmath at 1150 dps
(coefficients from stage1_exact_cz_nice.json, == M4a w4pp_closed_form.json term by term),
and compare against M4a_decomposition/W4pp_1060.txt (1060 digits, planner-checked).
Precision set AFTER all imports (dilog_exact/check_conversion traps do not apply here but be safe)."""
import sys, json, pickle, re, time
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/M4b2_w4pp_nice')
import m4b2_common as mc
from fractions import Fraction as F
import mpmath as mp
mp.mp.dps = 1150
mqf = lambda f: mp.mpf(f.numerator)/f.denominator

t0 = time.time()
d = json.load(open(str(mc.HERE/'stage1_exact_cz_nice.json')))
AT = pickle.load(open(str(mc.HERE/'stage2_atoms_nice.pkl'), 'rb'))['atoms']
cr = [F(x) for x in d['c_rat']]; cs = [F(x) for x in d['c_sq']]

def parse_qs2_mpf(s):
    m = re.fullmatch(r"\\((.+?)([+-])(\\d+(?:/\\d+)?)\\*sqrt2\\)", s)
    if m:
        q = F(m.group(3)); q = q if m.group(2) == '+' else -q
        return mqf(F(m.group(1))) + mp.sqrt(2)*mqf(q)
    if s.endswith('*sqrt2'): return mp.sqrt(2)*mqf(F(s[:-6]))
    return mqf(F(s))

def atom_value(enc):
    v = mp.mpf(1)
    for a in enc:
        if a[0] == 'log': v *= mp.log(parse_qs2_mpf(a[1]))
        elif a[0] == 'pi2': v *= mp.pi**2
        elif a[0] == 'z3': v *= mp.zeta(3)
        elif a[0] == 'li': v *= mp.polylog(int(a[1]), parse_qs2_mpf(a[2]))
    return v

Rpp = mp.mpf(0)
for k, (name, enc, _) in enumerate(AT):
    coef = mqf(cr[k]) + mp.sqrt(2)*mqf(cs[k])
    if coef:
        Rpp += coef*atom_value(enc)
W4 = mp.mpf(open(str(mc.THEORY/'M4a_decomposition/W4pp_1060.txt')).read().strip())
err = abs(Rpp - W4)
digits = -mp.log10(err) if err else mp.inf
mc.log(f"R''_nice (1150 dps) = {mp.nstr(Rpp, 80)}")
mc.log(f"W4'' (1060 digits)  = {mp.nstr(W4, 80)}")
mc.log(f"|R''_nice - W4''| = {mp.nstr(err, 3)}  -> agreement to {mp.nstr(digits, 6)} decimal digits")
assert digits > 1050, digits
# also: cross-check certificate coefficients == M4a w4pp_closed_form.json evaluated directly
m4a = json.load(open(str(mc.THEORY/'M4a_decomposition/w4pp_closed_form.json')))['closed_form_terms']
s = mp.mpf(0)
name2idx = {name: k for k, (name, _, _) in enumerate(AT)}
for coeff, expr in m4a:
    if expr.startswith('sqrt2*'): coef = mp.sqrt(2)*mqf(F(coeff)); name = expr[6:]
    else: coef = mqf(F(coeff)); name = expr
    s += coef*atom_value(AT[name2idx[name]][1])
err2 = abs(s - Rpp)
mc.log(f"cross-check: direct eval of M4a w4pp_closed_form.json - certificate R''_nice = {mp.nstr(err2, 3)}")
mc.log(f'step09 done s={time.time()-t0:.1f}')
'''
open(here + '/step09_sanity_1060.py', 'w').write(step09)
import subprocess
PY = base + '/.venv/bin/python'
r = subprocess.run(['nice', PY, 'step09_sanity_1060.py'], cwd=here, capture_output=True, text=True)
open(here + '/step09_sanity_1060.out', 'w').write(r.stdout + r.stderr)
print('rc=', r.returncode); print(r.stdout[-1200:])
if r.returncode: print('STDERR:', r.stderr[-1500:])

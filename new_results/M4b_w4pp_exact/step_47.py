import os, subprocess, time
Z3 = 'Zeta' + '(3)'
ROOT = os.path.join('<home>/Downloads/Research', Z3)
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
PY = os.path.join(ROOT, '.venv/bin/python')
step05 = r'''"""step05: numeric verification of R'' against W4'' (quad_90.json, 90 digits) and readable print.
R'' = sum_k (c_rat,k + sqrt2*c_sq,k) * B_k with c from stage1_exact_cz.json.
Evaluated directly with mpmath at 110 dps (precision set after all imports)."""
import sys, json, pickle, re, time
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/M4b_w4pp_exact')
import m4b_common as mc
from fractions import Fraction as F
import mpmath as mp
mp.mp.dps = 110

t0 = time.time()
d = json.load(open(str(mc.HERE/'stage1_exact_cz.json')))
AT = pickle.load(open(str(mc.HERE/'stage2_atoms.pkl'), 'rb'))['atoms']
cr = [F(x) for x in d['c_rat']]; cs = [F(x) for x in d['c_sq']]

def parse_qs2_mpf(s):
    m = re.fullmatch(r"\((.+?)([+-])(\d+(?:/\d+)?)\*sqrt2\)", s)
    if m:
        q = F(m.group(3)); q = q if m.group(2) == '+' else -q
        return mp.mpf(F(m.group(1))) + mp.sqrt(2)*mp.mpf(q)
    if s.endswith('*sqrt2'): return mp.sqrt(2)*mp.mpf(F(s[:-6]))
    return mp.mpf(s)

def atom_value(enc):
    v = mp.mpf(1)
    for a in enc:
        if a[0] == 'log': v *= mp.log(parse_qs2_mpf(a[1]))
        elif a[0] == 'pi2': v *= mp.pi**2
        elif a[0] == 'z3': v *= mp.zeta(3)
        elif a[0] == 'li': v *= mp.polylog(int(a[1]), parse_qs2_mpf(a[2]))
    return v

Rpp = mp.mpf(0); terms = []
for k, (name, enc, _) in enumerate(AT):
    coef = mp.mpf(cr[k]) + mp.sqrt(2)*mp.mpf(cs[k])
    if coef:
        Rpp += coef*atom_value(enc)
        terms.append((name, cr[k], cs[k]))
W4pp = mp.mpf(json.load(open(str(mc.THEORY/'M4a_decomposition/quad_90.json')))['W4pp'])
err = abs(Rpp - W4pp)
mc.log(f"R'' = {mp.nstr(Rpp, 100)}")
mc.log(f"W4pp= {mp.nstr(W4pp, 100)}")
mc.log(f"|R'' - W4''| = {mp.nstr(err, 3)}  (quad_90.json has 90 digits)")
assert err < mp.mpf('1e-89')

# readable print, Conjecture-4 style
def qs(c1, c2):
    """format c1 + c2*sqrt2"""
    if not c2: return f'{c1}'
    if not c1: return f'({c2})*sqrt2' if c2 != 1 else 'sqrt2'
    return f'({c1}{"+" if c2 > 0 else "-"}{abs(c2)}*sqrt2)'
li = [(n, a, b) for n, a, b in terms if n.startswith('Li4')]
pr = [(n, a, b) for n, a, b in terms if not n.startswith('Li4')]
lines = []
lines.append("R'' =")
for n, a, b in li:
    lines.append(f'    {"+" if (a, b) >= (0, 0) and (a > 0 or b > 0) else "-"}{""}')
out = []
for n, a, b in li + pr:
    coef = qs(abs(a), abs(b)) if (a, b) != (F(0), F(0)) else '0'
    sgn = '-' if (a < 0 or (a == 0 and b < 0)) else '+'
    out.append((sgn, coef, n))
msg = "R'' = " + ' '.join(f'{s} {c}*{n}' for s, c, n in out)
mc.log('readable form:')
for chunk in [out[i:i+4] for i in range(0, len(out), 4)]:
    mc.log('   ' + ' '.join(f'{s} {c}*{n}' for s, c, n in chunk))
mc.log(f'step05 done s={time.time()-t0:.1f}')
'''
open('step05_numeric_Rpp.py','w').write(step05)
r = subprocess.run([PY, 'step05_numeric_Rpp.py'], capture_output=True, text=True, cwd=WORK)
print('EXIT', r.returncode)
print(r.stdout)
print(r.stderr[-2000:] if r.returncode else '')

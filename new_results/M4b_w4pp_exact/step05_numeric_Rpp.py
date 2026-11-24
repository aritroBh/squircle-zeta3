"""step05: numeric verification of R'' against W4'' (quad_90.json, 90 digits) and readable print.
R'' = sum_k (c_rat,k + sqrt2*c_sq,k) * B_k with c from stage1_exact_cz.json.
Evaluated directly with mpmath at 110 dps (precision set after all imports)."""
import sys, json, pickle, re, time
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/M4b_w4pp_exact')
import m4b_common as mc
from fractions import Fraction as F
import mpmath as mp
mp.mp.dps = 110
mqf = lambda f: mp.mpf(f.numerator)/f.denominator

t0 = time.time()
d = json.load(open(str(mc.HERE/'stage1_exact_cz.json')))
AT = pickle.load(open(str(mc.HERE/'stage2_atoms.pkl'), 'rb'))['atoms']
cr = [F(x) for x in d['c_rat']]; cs = [F(x) for x in d['c_sq']]

def parse_qs2_mpf(s):
    m = re.fullmatch(r"\((.+?)([+-])(\d+(?:/\d+)?)\*sqrt2\)", s)
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

Rpp = mp.mpf(0); terms = []
for k, (name, enc, _) in enumerate(AT):
    coef = mqf(cr[k]) + mp.sqrt(2)*mqf(cs[k])
    if coef:
        Rpp += coef*atom_value(enc)
        terms.append((name, cr[k], cs[k]))
W4pp = mp.mpf(json.load(open(str(mc.THEORY/'M4a_decomposition/quad_90.json')))['W4pp'])
err = abs(Rpp - W4pp)
mc.log(f"R'' = {mp.nstr(Rpp, 100)}")
mc.log(f"W4pp= {mp.nstr(W4pp, 100)}")
mc.log(f"|R'' - W4''| = {mp.nstr(err, 3)}  (quad_90.json has 90 digits)")
assert err < mp.mpf('1e-89')

# readable print, Conjecture-4 style (coefficient (c_rat + c_sq*sqrt2) printed exactly)
def coef_str(a, b):
    if not b: return f'{a}'
    if not a: return f'{b}*sqrt2'
    return f'({a}{"+" if b > 0 else "-"}{abs(b)}*sqrt2)'
mc.log("readable form (R'' = sum (c_rat + c_sq*sqrt2) * atom):")
lines = [f'   {"+" if (a > 0 or (a == 0 and b >= 0)) else "-"} {coef_str(abs(a), abs(b))}*{n}'
         if False else f'   ({coef_str(a, b)}) * {n}' for n, a, b in terms]
for chunk in [lines[i:i+2] for i in range(0, len(lines), 2)]:
    mc.log(' ; '.join(chunk))
mc.log(f'step05 done s={time.time()-t0:.1f}')

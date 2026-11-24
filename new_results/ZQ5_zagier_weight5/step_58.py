import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json, math, pickle, random
import mpmath as mp
import sympy as sp
import w5lib as W
import pipe
D=5
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
vals = pipe.load_vals(D, 1150)
RM = json.load(open(pipe.cand_dir(D)+"/rel_min_atom.json"))
names = RM["names"]; c = [Fraction(x) for x in RM["c"]]
cmap = dict(zip(names, c))
N = len(cands); d = 3
prods = pipe.prod_values(C, vals)

# (b) exact symbol
T = {}
for i, cd in enumerate(cands):
    q = cmap.get(f"Li5({pipe.fmt_qd(cd['x'])})", Fraction(0))
    if not q: continue
    s = W.symbol_li5_dense(cd['vx'], cd['v1m'], d)
    for idx, cf in enumerate(s):
        if cf: T[idx] = T.get(idx, Fraction(0)) + q*cf
for (kind, msk, nm, pv) in prods:
    q = cmap.get(nm, Fraction(0))
    if not q or kind != "log5": continue
    v = W.prod_symbol_dense(msk, d)
    for idx, cf in enumerate(v):
        if cf: T[idx] = T.get(idx, Fraction(0)) + q*cf
Tz = {k:v for k,v in T.items() if v}
print("(b) exact full symbol of D=5 RHS:", "ZERO" if not Tz else f"NONZERO {Tz}")

# (c) cobrackets
def to_tup(idx):
    tup = []
    x = idx
    for _ in range(5): tup.append(x % d); x //= d
    return tuple(reversed(tup))
dic2 = {}
for idx, cf in T.items(): dic2[to_tup(idx)] = cf
for mode in ("full","cut"):
    r = W.cobracket32(dic2, d, mode)
    print(f"(c) delta_{{3,2}} total [{mode}]: {len(r)} nonzero (0 expected)")
# also Li5-part separately
Tl = {}
for i, cd in enumerate(cands):
    q = cmap.get(f"Li5({pipe.fmt_qd(cd['x'])})", Fraction(0))
    if not q: continue
    s = W.symbol_li5_dense(cd['vx'], cd['v1m'], d)
    for idx, cf in enumerate(s):
        if cf: Tl[idx] = Tl.get(idx, Fraction(0)) + q*cf
dicl = {}
for idx, cf in Tl.items(): dicl[to_tup(idx)] = cf
for mode in ("full","cut"):
    r = W.cobracket32(dicl, d, mode)
    print(f"(c) delta_{{3,2}} Li5-part [{mode}]: {len(r)} nonzero (0 expected)")

# (d) factorisations
print("(d) used Li5 arguments:")
for i, cd in enumerate(cands):
    q = cmap.get(f"Li5({pipe.fmt_qd(cd['x'])})", Fraction(0))
    if not q: continue
    print(f"  x={pipe.fmt_qd(cd['x']):>18s} coeff={str(q):>10s}  v(x)={cd['vx']}  v(1-x)={cd['v1m']}")

# (e)
used = [(nm, cf) for nm, cf in cmap.items() if nm != "T" and cf]
dens = [x.denominator for _,x in used]; nums = [abs(x.numerator) for _,x in used]
print("(e) terms:", len(used), " max|num| =", max(nums), " lcm den =", math.lcm(*dens), "=", sp.factorint(math.lcm(*dens)))

# fresh independent evaluation at 1400 digits
P = 1400
mp.mp.dps = P + 30
lphi = mp.log((1+mp.sqrt(5))/2); ls5 = mp.log(mp.sqrt(5)); l2 = mp.log(2); pi2 = mp.pi**2
tot = W.Lchi_hurwitz(D, P)/mp.sqrt(5)
for nm, cf in cmap.items():
    if nm == "T" or not cf: continue
    if nm == "z5": v = mp.zeta(5)
    elif nm.startswith("Li5("):
        i = [j for j, cd in enumerate(cands) if f"Li5({pipe.fmt_qd(cd['x'])})" == nm][0]
        v = W.li5_series(cands[i]['x'].emb(P+30), P)
    else:
        v = mp.mpf(1)
        for f in nm.split("*"): v *= {"lphi": lphi, "ls5": ls5, "l2": l2, "pi2": pi2}[f]
    tot -= cf*v
print("fresh 1400-digit residual T - RHS:", mp.nstr(tot, 4))
'''
open('steps/step_33_d5_verify.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_33_d5_verify.py'],capture_output=True,text=True,timeout=2400)
out = r.stdout
lines = [l for l in out.splitlines() if not l.startswith("D=5 LLL")]
print("\n".join(lines)); print("STDERR:", r.stderr[-1000:])
with open('progress.log','a') as f: f.write("step_33_d5_verify: (b)-(e) + fresh 1400d\n")

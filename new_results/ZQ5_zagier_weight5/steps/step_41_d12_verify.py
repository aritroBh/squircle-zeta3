
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json, math, pickle
import mpmath as mp
import sympy as sp
import w5lib as W
import pipe
D=12
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
vals = pipe.load_vals(D, 1150)
RM = json.load(open(pipe.cand_dir(D)+"/rel_min_atom.json"))
names = RM["names"]; c = [Fraction(x) for x in RM["c"]]
cmap = dict(zip(names, c))
N = len(cands); d = 3
prods = pipe.prod_values(C, vals)

# (a) two-precision span check using saved atom relations is already done (identical: False);
# redo at kernel level with span test:
K = [list(map(int, r)) for r in json.load(open(pipe.cand_dir(D)+"/kernel.json"))["K"]]
R8 = pipe.phase_lll(D, C, cands, K, vals, P=830, use_atoms=False)
R5 = pipe.phase_lll(D, C, cands, K, vals, P=520, use_atoms=False)
a8 = R8["exact_rels"][0][1]; a5 = R5["exact_rels"][0][1]
Zl = [zz for _,zz,_ in R8["exact_zeros"]]
Zm = sp.Matrix(Zl)
print("(a) two-precision: a0s", a5[0], a8[0], " (a5-a8) in zero span:",
      sp.Matrix(Zl+[[a5[i]-a8[i] for i in range(len(a5))]]).rank() == Zm.rank())

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
print("(b) exact full symbol of D=12 RHS:", "ZERO" if not Tz else f"NONZERO {len(Tz)} entries")

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
    print(f"(c) delta_{{3,2}} total [{mode}]: {len(r)} nonzero")

# (d) factorisations (compact: count + print all)
used = [(nm, cf) for nm, cf in cmap.items() if nm != "T" and cf]
li5used = [(nm, cf) for nm, cf in used if nm.startswith("Li5")]
print("(d) Li5 terms:", len(li5used), "; every used x has x and 1-x S-units with coords:")
for nm, cf in li5used:
    i = [j for j, cd in enumerate(cands) if f"Li5({pipe.fmt_qd(cd['x'])})" == nm][0]
    cd = cands[i]
    print(f"  x={pipe.fmt_qd(cd['x']):>18s} coeff={str(cf):>11s}  v(x)={cd['vx']}  v(1-x)={cd['v1m']}")

# (e)
dens = [x.denominator for _,x in used]; nums = [abs(x.numerator) for _,x in used]
print("(e) terms:", len(used), " max|num| =", max(nums), " lcm den =", math.lcm(*dens), "=", sp.factorint(math.lcm(*dens)))

# fresh 1400-digit residual
P = 1400
mp.mp.dps = P + 30
lE3 = mp.log(2+mp.sqrt(3)); ls3 = mp.log(mp.sqrt(3)); lu = mp.log(1+mp.sqrt(3)); pi2 = mp.pi**2
tot = W.Lchi_hurwitz(D, P)/mp.sqrt(3)
for nm, cf in cmap.items():
    if nm == "T" or not cf: continue
    if nm == "z5": v = mp.zeta(5)
    elif nm.startswith("Li5("):
        i = [j for j, cd in enumerate(cands) if f"Li5({pipe.fmt_qd(cd['x'])})" == nm][0]
        v = W.li5_series(cands[i]['x'].emb(P+30), P)
    else:
        v = mp.mpf(1)
        for f in nm.split("*"): v *= {"lE3": lE3, "ls3": ls3, "lu": lu, "pi2": pi2}[f]
    tot -= cf*v
print("fresh 1400-digit residual T - RHS:", mp.nstr(tot, 4))

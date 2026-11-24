
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json, math, time
import mpmath as mp
import sympy as sp
import w5lib as W
import pipe

D=8
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
K = [list(map(int, r)) for r in json.load(open(pipe.cand_dir(D)+"/kernel.json"))["K"]]
vals = pipe.load_vals(D, 1150)
RM = json.load(open(pipe.cand_dir(D)+"/rel_min.json"))
names = RM["names"]; c = [Fraction(x) for x in RM["c"]]
cmap = dict(zip(names, c))
N = len(cands)
d = 2

# (a) two-precision: exact relation vectors identical
R5 = pipe.phase_lll(D, C, cands, K, vals, P=520, use_atoms=True)
R8 = pipe.phase_lll(D, C, cands, K, vals, P=830, use_atoms=True)
a5 = R5["exact_rels"][0][1]; a8 = R8["exact_rels"][0][1]
print("(a) exact rel at 520:", a5[:6], "... resid", mp.nstr(R5['exact_rels'][0][2],2))
print("    exact rel at 830:", a8[:6], "... resid", mp.nstr(R8['exact_rels'][0][2],2))
print("    identical vectors:", a5 == a8)
# canonical representatives: expand each on atoms, then check difference in zero span
diff58 = [a5[i]-a8[i] for i in range(len(a5))]
Zl5 = [zz for _, zz, _ in R8["exact_zeros"]]
Zm5 = sp.Matrix([[Fraction(x) for x in zz] for zz in Zl5])
print("    (a5-a8) in Q-span of exact zeros:", sp.Matrix(Zl5 + [diff58]).rank() == Zm5.rank())

# (b) exact full symbol of the minimised RHS
T = {}
for i, cd in enumerate(cands):
    q = cmap.get(f"Li5({pipe.fmt_qd(cd['x'])})", Fraction(0))
    if not q: continue
    s = W.symbol_li5_dense(cd['vx'], cd['v1m'], d)
    for idx, cf in enumerate(s):
        if cf: T[idx] = T.get(idx, Fraction(0)) + q*cf
prods = pipe.prod_values(C, vals)
Tprod = {}
for (kind, msk, nm, pv) in prods:
    q = cmap.get(nm, Fraction(0))
    if not q or kind != "log5": continue
    v = W.prod_symbol_dense(msk, d)
    for idx, cf in enumerate(v):
        if cf: T[idx] = T.get(idx, Fraction(0)) + q*cf
        if cf: Tprod[idx] = Tprod.get(idx, Fraction(0)) + q*cf
Tli5 = {}
for i, cd in enumerate(cands):
    q = cmap.get(f"Li5({pipe.fmt_qd(cd['x'])})", Fraction(0))
    if not q: continue
    s = W.symbol_li5_dense(cd['vx'], cd['v1m'], d)
    for idx, cf in enumerate(s):
        if cf: Tli5[idx] = Tli5.get(idx, Fraction(0)) + q*cf
Tz = {k:v for k,v in T.items() if v}
print("(b) exact full symbol of RHS (Li5 + log^5 products):", "ZERO" if not Tz else f"NONZERO {Tz}")

# (c) cobracket of total, Li5-part, products-part, both modes
def sparse_dic(dense):
    dic = {}
    for idx, cf in enumerate(dense):
        if cf:
            tup = []
            x = idx
            for _ in range(5): tup.append(x % d); x //= d
            dic[tuple(reversed(tup))] = dic.get(tuple(reversed(tup)), Fraction(0)) + cf
    return dic
for tag, TT in [("total", T), ("Li5-part", Tli5), ("products-part", Tprod)]:
    dic = {k: v for k, v in TT.items()}
    # need tuple keys: rebuild
    dic2 = {}
    for idx, cf in TT.items():
        tup = []
        x = idx
        for _ in range(5): tup.append(x % d); x //= d
        dic2[tuple(reversed(tup))] = cf
    r_full = W.cobracket32(dic2, d, "full")
    r_cut = W.cobracket32(dic2, d, "cut")
    print(f"(c) delta_{{3,2}} [{tag}]: full={len(r_full)} nonzero, cut={len(r_cut)} nonzero (0 expected)")

# (d) factorisations of used Li5 args
print("(d) used Li5 arguments: x and 1-x as S-units (letter coords", C["names"], "):")
for i, cd in enumerate(cands):
    q = cmap.get(f"Li5({pipe.fmt_qd(cd['x'])})", Fraction(0))
    if not q: continue
    print(f"  x={pipe.fmt_qd(cd['x']):>18s} coeff={str(q):>10s}  v(x)={cd['vx']}  v(1-x)={cd['v1m']}")

# (e) lcm / max numerator
used = [(nm, cf) for nm, cf in cmap.items() if nm != "T" and cf]
dens = [x.denominator for _, x in used]; nums = [abs(x.numerator) for _, x in used]
print("(e) terms:", len(used), " max|num| =", max(nums), " lcm den =", math.lcm(*dens),
      "= factor:", sp.factorint(math.lcm(*dens)))

# comparison with M4i: (theirs - mine) in zero span?
M4 = json.load(open("../M4i_w5_noL/L_relation.json"))["terms"]
m4c = {}
for key, (p, q) in M4.items():
    m4c[key] = Fraction(q)
# express their vector on MY names (products converted lam->2*ls2)
theirs = {}
for key, cf in m4c.items():
    if key == "z5": theirs["z5"] = cf
    elif key.startswith("Li5("):
        inner = key[4:-1]
        vv = eval(inner, {"sqrt2": mp.sqrt(2)})
        bi = min(range(N), key=lambda i: abs(cands[i]['x'].emb(150)-vv))
        theirs[f"Li5({pipe.fmt_qd(cands[bi]['x'])})"] = cf
    else:
        fs = key.split("*")
        na = fs.count("lam"); nb = fs.count("L"); npi = fs.count("pi2")
        key2 = "*".join(["ls2"]*na + ["lE"]*nb + ["pi2"]*npi)
        theirs[key2] = cf * 2**na
Zl = [zz for _, zz, _ in R8["exact_zeros"]]
Zf = [[Fraction(x) for x in zz] for zz in Zl]
diff = []
for nm in names:
    diff.append(theirs.get(nm, Fraction(0)) - cmap.get(nm, Fraction(0)))
Zm = sp.Matrix(Zf)
print("M4i - mine in Q-span of exact zeros:", sp.Matrix(Zf + [diff]).rank() == Zm.rank())

# fresh independent evaluation at 1200 digits
P = 1200
mp.mp.dps = P + 30
lam = mp.log(mp.sqrt(2)); lL = mp.log(1+mp.sqrt(2)); pi2 = mp.pi**2
tot = W.Lchi_hurwitz(D, P)/mp.sqrt(2)   # T
for nm, cf in cmap.items():
    if nm == "T" or not cf: continue
    if nm == "z5": v = mp.zeta(5)
    elif nm.startswith("Li5("):
        i = [j for j, cd in enumerate(cands) if f"Li5({pipe.fmt_qd(cd['x'])})" == nm][0]
        v = W.li5_series(cands[i]['x'].emb(P+30), P)
    else:
        v = mp.mpf(1)
        for f in nm.split("*"): v *= {"ls2": lam, "lE": lL, "pi2": pi2}[f]
    tot -= cf*v
print("fresh 1200-digit residual T - RHS:", mp.nstr(tot, 4))

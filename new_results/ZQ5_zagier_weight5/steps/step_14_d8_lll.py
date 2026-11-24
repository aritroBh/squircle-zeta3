
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json, time, math
import mpmath as mp
import sympy as sp
import flint
import w5lib as W
import pipe

D = 8
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
K = [list(map(int, r)) for r in json.load(open(pipe.cand_dir(D)+"/kernel.json"))["K"]]
vals = pipe.phase_numerics(D, C, cands, P=1150)

def get_rel(P, use_atoms):
    return pipe.phase_lll(D, C, cands, K, vals, P=P, use_atoms=use_atoms)

t0=time.time()
R1a = get_rel(520, True)
R2a = get_rel(830, True)
R1k = get_rel(520, False)
R2k = get_rel(830, False)
print(f"total LLL time {time.time()-t0:.0f}s")

def prim(a):
    g = 0
    for x in a: g = math.gcd(g, abs(x))
    a = [x//g for x in a]
    if a[0] < 0: a = [-x for x in a]
    return a

for tag, R in [("atoms", (R1a, R2a)), ("kernel", (R1k, R2k))]:
    r1, r2 = R
    b1 = r1["rels"][0][1]; b2 = r2["rels"][0][1]
    m1 = r1["rels"][0][0]; m2 = r2["rels"][0][0]
    res1 = r1["rels"][0][2]; res2 = r2["rels"][0][2]
    p1 = prim(b1); p2 = prim(b2)
    same = p1 == p2
    print(tag, ": best rel max|a| P520=", m1, " P830=", m2, " identical after primitive: ", same)
    print("   resid P520=", mp.nstr(res1,2), " resid P830=", mp.nstr(res2,2))
    if not same:
        print("   P520:", p1)
        print("   P830:", p2)

# use atom-level 830 relation as main
names = R2a["names"]; a = R2a["rels"][0][1]
a0 = a[0]
coef = [Fraction(-x, a0) for x in a]   # T = sum coef[i]*v_i, coef[0]=1
assert coef[0] == 1
# coefficient vector on Li5 atoms + products (skip T, z5 handled separately)
li5_coef = coef[2:2+len(cands)]
z5_coef = coef[1]
prod_coef = coef[2+len(cands):]

# compare with M4i
M4 = json.load(open("../M4i_w5_noL/L_relation.json"))["terms"]
m4vec = {}
for key, (p, q) in M4.items():
    m4vec[key] = Fraction(q)
# build my vector in M4i naming
def mykey_li5(c): return f"Li5({pipe.fmt_qd(c['x'])})"
mine = {}
for c, cf in zip(cands, li5_coef):
    if cf: mine[mykey_li5(c)] = cf
mine["z5"] = z5_coef
prods = R2a["prods"]
for (kind, msk, nm), cf in zip(prods, prod_coef):
    if not cf: continue
    # convert my letters to M4i letters: ls2 -> lam/2  (coeff *= 2^{-#ls2})
    if kind == "log5":
        na = msk.count(0); nb = msk.count(1)
        key = "*".join(["lam"]*na + ["L"]*nb)
        cf = cf / 2**na
    elif kind == "pi2log3":
        na = msk.count(0); nb = msk.count(1)
        key = "*".join(["lam"]*na + ["L"]*nb + ["pi2"])
        cf = cf / 2**na
    else:
        i = msk[0]
        key = "*".join((["lam"] if i==0 else ["L"]) + ["pi2","pi2"])
        cf = cf / (2 if i==0 else 1)
    mine[key] = mine.get(key, Fraction(0)) + cf
allk = sorted(set(mine) | set(m4vec))
print("--- comparison with M4i (mine vs theirs) ---")
diff_any = False
for k in allk:
    m1 = mine.get(k, Fraction(0)); m2 = m4vec.get(k, Fraction(0))
    flag = "" if m1 == m2 else "   <<< DIFF"
    if m1 != m2: diff_any = True
    print(f"{k:>26s}  {str(m1):>14s}  {str(m2):>14s}{flag}")
print("IDENTICAL to M4i:", not diff_any)
json.dump(dict(a=a, names=names), open(pipe.cand_dir(D)+"/rel_raw.json","w"))

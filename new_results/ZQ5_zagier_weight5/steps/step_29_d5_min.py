
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json, math, time, itertools
import mpmath as mp
import sympy as sp
import w5lib as W
import pipe
D=5
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
K = [list(map(int, r)) for r in json.load(open(pipe.cand_dir(D)+"/kernel.json"))["K"]]
vals = pipe.load_vals(D, 1150)
N = len(cands)

R8 = pipe.phase_lll(D, C, cands, K, vals, P=830, use_atoms=False)
R5 = pipe.phase_lll(D, C, cands, K, vals, P=520, use_atoms=False)
a8 = R8["exact_rels"][0][1]; a5 = R5["exact_rels"][0][1]
print("(a) identical exact relation at 520 and 830:", a5 == a8)
if a5 != a8:
    Zl = [zz for _,zz,_ in R8["exact_zeros"]]
    Zm = sp.Matrix(Zl)
    print("    (a5-a8) in zero span:", sp.Matrix(Zl+[[a5[i]-a8[i] for i in range(len(a5))]]).rank() == Zm.rank())

# rebuild keep/dropped mapping (same logic as phase_lll)
mp.mp.dps = 1180
li5v = [vals[f"li5:{i}"] for i in range(N)]
combo_vals = [sum(mp.mpf(c)*v for c,v in zip(row, li5v)) for row in K]
keep_idx = [j for j,v in enumerate(combo_vals) if abs(v) >= mp.mpf(10)**(-(830-80))]
drop_idx = [j for j,v in enumerate(combo_vals) if abs(v) < mp.mpf(10)**(-(830-80))]
names = R8["names"]; prods = R8["prods"]
nk = len(keep_idx)
a0 = a8[0]
print("a0 =", a0, "=", sp.factorint(abs(a0)), " kept combos:", nk, " dropped:", drop_idx)

# atom-space coefficient vector on [z5] + Li5(N) + products(M)
M = len(prods)
c = [Fraction(-x, a0) for x in a8]
c_z5 = c[1]
c_li5 = [Fraction(0)]*N
for j in range(nk):
    if c[2+j]:
        Kj = K[keep_idx[j]]
        for k in range(N):
            if Kj[k]: c_li5[k] += c[2+j]*Kj[k]
c_prod = c[2+nk:]
cvec = [c_z5] + c_li5 + c_prod

# zero lattice on atom space
Z = []
for _, zz, _ in R8["exact_zeros"]:
    zv = [Fraction(zz[1])] + [Fraction(0)]*N + [Fraction(x) for x in zz[2+nk:]]
    for j in range(nk):
        if zz[2+j]:
            Kj = K[keep_idx[j]]
            for k in range(N):
                if Kj[k]: zv[1+k] += Fraction(zz[2+j])*Kj[k]
    Z.append(zv)
for j in drop_idx:
    zv = [Fraction(0)] + [Fraction(x) for x in K[j]] + [Fraction(0)]*M
    Z.append(zv)
Zm = sp.Matrix(Z)
print("atom zero lattice rows:", len(Z), " rank:", Zm.rank())

def cost(cv):
    return max(abs(x.numerator) for x in cv)

# greedy: eliminate worst Li5 coordinate using best zero vector
cv = list(cvec)
Zlive = [list(z) for z in Z]
t0=time.time()
for it in range(400):
    cur = cost(cv)
    best = None
    li5span = range(1, 1+N)
    for i in li5span:
        if cv[i] == 0: continue
        for zi, z in enumerate(Zlive):
            if z[i] == 0: continue
            t = -cv[i]/z[i]
            cand = [cv[k]+t*z[k] for k in range(len(cv))]
            co = cost(cand)
            if best is None or co < best[0]:
                best = (co, i, zi, t, cand)
    if best and best[0] < cur:
        cv = best[4]
    else:
        break
print(f"greedy done ({time.time()-t0:.0f}s) cost {cost(cv)}")
used = []
nm_full = ["z5"] + [f"Li5({pipe.fmt_qd(cd['x'])})" for cd in cands] + [p[2] for p in prods]
for nm, cf in zip(nm_full, cv):
    if cf: used.append((nm, cf))
dens = [x.denominator for _,x in used]; nums = [abs(x.numerator) for _,x in used]
print("terms:", len(used), " max|num| =", max(nums), " lcm den =", math.lcm(*dens), "=", sp.factorint(math.lcm(*dens)))
for nm, cf in used: print(f"  {nm:>26s}  {cf}")
json.dump(dict(names=["T"]+nm_full, c=[str(Fraction(-1) if i==0 else x) for i,x in enumerate([Fraction(-1)]+list(cv))]),
          open(pipe.cand_dir(D)+"/rel_min.json","w"))

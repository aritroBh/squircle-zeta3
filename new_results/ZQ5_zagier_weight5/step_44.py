import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json
import mpmath as mp
import sympy as sp
import w5lib as W
import pipe
D=8
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
vals = pipe.load_vals(D, 1150)
K = [list(map(int, r)) for r in json.load(open(pipe.cand_dir(D)+"/kernel.json"))["K"]]
N = len(cands); d = 2
prods = pipe.prod_values(C, vals)

def symbol_of_vec(cmap):
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
    return {k: v for k, v in T.items() if v}

# 1. M4i's RHS symbol with my code
M4 = json.load(open("../M4i_w5_noL/L_relation.json"))["terms"]
theirs = {}
for key, (p, q) in M4.items():
    cf = Fraction(q)
    if key == "z5": theirs["z5"] = cf
    elif key.startswith("Li5("):
        vv = eval(key[4:-1], {"sqrt2": mp.sqrt(2)})
        bi = min(range(N), key=lambda i: abs(cands[i]['x'].emb(150)-vv))
        theirs[f"Li5({pipe.fmt_qd(cands[bi]['x'])})"] = cf
    else:
        fs = key.split("*")
        na = fs.count("lam"); nb = fs.count("L"); npi = fs.count("pi2")
        key2 = "*".join(["ls2"]*na + ["lE"]*nb + ["pi2"]*npi)
        theirs[key2] = cf * 2**na
S1 = symbol_of_vec(theirs)
print("1. symbol of M4i RHS with MY code:", "ZERO" if not S1 else f"NONZERO ({len(S1)} entries) e.g. {dict(list(S1.items())[:4])}")

# 2. symbols of the 5 exact-zero rows
R8 = pipe.phase_lll(D, C, cands, K, vals, P=830, use_atoms=True)
names = R8["names"]
for j,(mx,zz,res) in enumerate(R8["exact_zeros"]):
    zmap = {names[i]: Fraction(x) for i,x in enumerate(zz) if x and names[i]!="T"}
    Sz = symbol_of_vec(zmap)
    print(f"2. symbol of exact-zero #{j} (max|a|={mx}):", "ZERO" if not Sz else f"NONZERO ({len(Sz)}) e.g. {dict(list(Sz.items())[:3])}")

# 3. value of M4i-RHS - my-RHS at 1200 digits
RM = json.load(open(pipe.cand_dir(D)+"/rel_min.json"))
cmap = dict(zip(RM["names"], [Fraction(x) for x in RM["c"]]))
mp.mp.dps = 1250
lam = mp.log(mp.sqrt(2)); lL = mp.log(1+mp.sqrt(2)); pi2 = mp.pi**2
li5f = {}
for i, cd in enumerate(cands):
    li5f[f"Li5({pipe.fmt_qd(cd['x'])})"] = W.li5_series(cd['x'].emb(1250), 1200)
def evalmap(cm):
    tot = mp.mpf(0)
    for nm, cf in cm.items():
        if nm == "T" or not cf: continue
        if nm == "z5": v = mp.zeta(5)
        elif nm.startswith("Li5("): v = li5f[nm]
        else:
            v = mp.mpf(1)
            for f in nm.split("*"): v *= {"ls2": lam, "lE": lL, "pi2": pi2}[f]
        tot += cf*v
    return tot
vdiff = evalmap(theirs) - evalmap(cmap)
print("3. value(M4i RHS - my RHS) at 1200d:", mp.nstr(vdiff, 3))
# 4. and M4i-mine coefficient diff vs zero span, using the EXPANDED zero lattice:
#    add (M4i-mine) as candidate new zero vector
'''
open('steps/step_22_symtest.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_22_symtest.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
with open('progress.log','a') as f: f.write("step_22_symtest: M4i symbol + zero-row symbols + value diff\n")

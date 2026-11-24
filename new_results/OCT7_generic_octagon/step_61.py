code = r'''
import mpmath as mp, sympy as sp, pickle, time, itertools
from fractions import Fraction
mp.mp.dps = 300
t_start = time.time()

# reload basis values (exact rational points) at 300 digits
l2, l3, l5 = mp.log(2), mp.log(3), mp.log(5)
z2 = mp.pi**2/6; z3 = mp.zeta(3)
logs = [(l2,'l2'),(l3,'l3'),(l5,'l5')]
def orbit(x): return {x, 1-x, 1/x, x/(x-1), 1/(1-x), (x-1)/x}
x1 = Fraction(-1,3); x2 = Fraction(-1,9)
seeds = [x1, x2, Fraction(-13,27), Fraction(-1,5), Fraction(1,6), Fraction(-7,9), Fraction(-19,81),
         Fraction(-1,2), Fraction(-2,3), Fraction(1,2)]
pts = set()
for s in seeds: pts |= orbit(s)
pts = sorted([p for p in pts if p < 1 and p != 0])
def f2m(q): return mp.mpf(q.numerator)/q.denominator
basis = [(mp.mpf(1),'1')] + logs + [(z2,'z2')]
for p in pts: basis.append((mp.polylog(2, f2m(p)), f'Li2({p})'))
for (a,na),(b,nb) in itertools.combinations_with_replacement(logs, 2): basis.append((a*b, na+nb))
basis.append((z3,'z3'))
for p in pts: basis.append((mp.polylog(3, f2m(p)), f'Li3({p})'))
for a,na in logs: basis.append((z2*a, 'z2'+na))
for p in pts:
    for a,na in logs: basis.append((mp.polylog(2, f2m(p))*a, f'Li2({p}){na}'))
for (a,na),(b,nb),(c,nc) in itertools.combinations_with_replacement(logs, 3): basis.append((a*b*c, na+nb+nc))
n0 = len(basis); vals = [v for v,_ in basis]; names = [nm for _,nm in basis]
print("full basis size:", n0)

# ---- quotient internal relations: LLL, collect candidate relation rows, verify at 700, row-reduce ----
D = sp.Integer(10)**280
def lattice_lll(vals):
    n = len(vals)
    rows = []
    for i in range(n):
        r = [sp.Integer(0)]*n + [sp.Integer(int(mp.floor(vals[i]*D + mp.mpf('0.5'))))]
        r[i] = sp.Integer(1)
        rows.append(r)
    return sp.Matrix(rows).lll(delta=sp.Rational(3,4)).tolist()

t0=time.time()
Lr = lattice_lll(vals)
print(f"[{time.time()-t0:.0f}s] LLL done")
cands = []
for r in Lr:
    if all(x==0 for x in r[:-1]): continue
    if abs(int(r[-1])) < 10**230:
        res = mp.fsum([int(r[i])*vals[i] for i in range(len(vals))])
        if abs(res) < mp.mpf(10)**(-260):
            cands.append(r)
print("candidate internal relations:", len(cands))
# verify each at 700 digits
mp.mp.dps = 700
l2h, l3h, l5h = mp.log(2), mp.log(3), mp.log(5)
z2h, z3h = mp.pi**2/6, mp.zeta(3)
lm = {'l2':l2h,'l3':l3h,'l5':l5h}
import re as _re
def bval(name):
    if name=='1': return mp.mpf(1)
    if name in lm: return lm[name]
    if name=='z2': return z2h
    if name=='z3': return z3h
    val = mp.mpf(1); rest = name
    if rest.startswith('z2'): val*=z2h; rest=rest[2:]
    mm = _re.match(r'^Li([23])\((-?\d+/\d+|-?\d+)\)', rest)
    if mm:
        q = Fraction(mm.group(2)); val *= mp.polylog(int(mm.group(1)), f2m(q)); rest = rest[mm.end():]
    for ch in _re.findall(r'l[235]', rest): val *= lm[ch]
    return val
vals700 = [bval(nm) for nm in names]
gen_rels = []
for r in cands:
    res = mp.fsum([int(r[i])*vals700[i] for i in range(len(vals700))])
    if abs(res) < mp.mpf(10)**(-670):
        gen_rels.append(r)
print("genuine internal relations (verified 700 digits):", len(gen_rels))
mp.mp.dps = 300
# row-reduce over Q to find pivots
R = sp.Matrix([[sp.Integer(int(x)) for x in r[:-1]] for r in gen_rels]) if gen_rels else sp.Matrix(0, len(vals), [])
if gen_rels:
    Rr, pivots = R.rref()
    pivots = list(pivots)
else:
    pivots = []
print("independent relations:", len(pivots), "-> drop", len(pivots), "basis elements")
keep = [i for i in range(len(vals)) if i not in set(pivots)]
indep = [(vals[i], names[i]) for i in keep]
print("independent basis size:", len(indep))
pickle.dump({'indep_names':[nm for _,nm in indep], 'gen_rels': [[int(x) for x in r[:-1]] for r in gen_rels],
             'names': names, 'pivots': pivots}, open('basis_indep.pkl','wb'))
with open("progress.log","a") as f:
    f.write(f"STEP 07: quotiented basis {n0} -> {len(indep)} independent ({len(gen_rels)} classical relations, verified 700d)\n")
print(f"total {time.time()-t_start:.0f}s")
'''
with open("step_07.py","w") as f: f.write(code)
import subprocess, sys
r = subprocess.run([sys.executable, "step_07.py"], capture_output=True, text=True, cwd="research/08_theory/OCT7_generic_octagon", timeout=3600)
print(r.stdout)
if r.stderr: print("STDERR:", r.stderr[-800:])

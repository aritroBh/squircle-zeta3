import os, sys, time, pickle
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
import mpmath as mp
mp.mp.dps = 250
from fractions import Fraction as F

r2 = mp.sqrt(2); X = mp.log(2); Y = mp.log(1+r2); b = r2-1
def q2f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r2
av = pickle.load(open("atom_vals.pkl","rb"))
Li5v = {k: mp.mpf(v) for k,v in av['Li5v'].items()}
Li4v = {k: mp.mpf(v) for k,v in av['Li4v'].items()}
Li3v = {k: mp.mpf(v) for k,v in av['Li3v'].items()}
Li2v = {k: mp.mpf(v) for k,v in av['Li2v'].items()}
L32_b1 = mp.mpf(av['L32_b1']); L32_h1 = mp.mpf(av['L32_h1'])
z5 = mp.zeta(5); z3 = mp.zeta(3); z2 = mp.pi**2/6

# map solution terms to values
data = pickle.load(open("symbol_fit2.pkl","rb"))
cols = data['cols']
xsol = pickle.load(open("symbol_fit2_solution.pkl","rb"))['x']
nz = [(label, kind, xv) for (label, kind), xv in zip(cols, xsol) if xv.p != 0 or xv.q != 0]

import re
def term_value(label, kind):
    m = re.match(r"Li5\((.*)\)", label)
    if m: return Li5v[m.group(1)]
    m = re.match(r"Li4\((.*)\)\*([XY])", label)
    if m: return Li4v[m.group(1)] * (X if m.group(2)=='X' else Y)
    m = re.match(r"Li3\((.*)\)\*X\^2/2", label)
    if m: return Li3v[m.group(1)] * X**2/2
    m = re.match(r"Li3\((.*)\)\*X\*Y", label)
    if m: return Li3v[m.group(1)] * X*Y
    m = re.match(r"Li3\((.*)\)\*Y\^2/2", label)
    if m: return Li3v[m.group(1)] * Y**2/2
    m = re.match(r"Li2\((.*)\)\*X\^3/6", label)
    if m: return Li2v[m.group(1)] * X**3/6
    m = re.match(r"Li2\((.*)\)\*X\^2Y/2", label)
    if m: return Li2v[m.group(1)] * X**2*Y/2
    m = re.match(r"Li2\((.*)\)\*XY\^2/2", label)
    if m: return Li2v[m.group(1)] * X*Y**2/2
    m = re.match(r"Li2\((.*)\)\*Y\^3/6", label)
    if m: return Li2v[m.group(1)] * Y**3/6
    m = re.match(r"Li3\((.*)\)\*Li2\((.*)\)", label)
    if m: return Li3v[m.group(1)] * Li2v[m.group(2)]
    m = re.match(r"Li2\((.*)\)\*Li2\((.*)\)\*([XY])", label)
    if m: return Li2v[m.group(1)] * Li2v[m.group(2)] * (X if m.group(3)=='X' else Y)
    m = re.match(r"Li32\(b,1\)", label)
    if m: return L32_b1
    m = re.match(r"Li32\(1/sqrt2,1\)", label)
    if m: return L32_h1
    raise ValueError(label)

W5 = mp.mpf(pickle.load(open("w5_250.pkl","rb"))['W5'])
cand = mp.mpf(0)
for label, kind, xv in nz:
    cand += q2f(xv) * term_value(label, kind)
R = W5 - cand
print("residual R =", mp.nstr(R, 30))
print("|R| =", mp.nstr(abs(R), 4))
# R should be symbol-invisible: fit over invisible basis
invis = {'z5': z5, 'z2*z3': z2*z3, 'z3*X^2': z3*X**2, 'z3*X*Y': z3*X*Y, 'z3*Y^2': z3*Y**2,
         'z2*X^3': z2*X**3, 'z2*X^2Y': z2*X**2*Y, 'z2*XY^2': z2*X*Y**2, 'z2*Y^3': z2*Y**3,
         'z2^2*X': z2**2*X, 'z2^2*Y': z2**2*Y}
# LLL PSLQ over Q(sqrt2): R = sum gamma_i * atom_i, gamma in Q(sqrt2)
from flint import fmpz_mat
names = list(invis)
vec = [R] + [invis[n] for n in names]
big = []
for v in vec: big += [v, r2*v]
def pslq_q2(target_and_basis, dps):
    n = len(target_and_basis)*2
    C = mp.mpf(10)**int(dps*0.9)
    rows = []
    for i in range(n):
        row = [0]*n + [int(mp.nint(C*target_and_basis[i//2][i%2]))] if False else None
    return None
# build lattice properly: unknowns: target coeff (must be nonzero), then for each basis atom a p,q pair
vals = []
for v in vec:
    vals += [v, r2*v]
n = len(vals)
C = mp.mpf(10)**220
rows = []
for i in range(n):
    row = [0]*(n+1)
    row[i] = 1
    row[n] = int(mp.nint(C*vals[i]))
    rows.append(row)
M = fmpz_mat(rows).lll()
cands = []
for i in range(n):
    v = [int(M[i,j]) for j in range(n)]
    if v[0] == 0: continue
    res = abs(sum(mp.mpf(v[j])*vals[j] for j in range(n)))
    maxc = max(abs(xx) for xx in v)
    cands.append((maxc, res, v))
cands.sort(key=lambda z: z[0])
for maxc, res, v in cands[:3]:
    print("maxc=%.2e res=%s" % (maxc, mp.nstr(res,3)))
    if res < mp.mpf('1e-60'):
        out = {}
        v0 = v[0]
        for i, nm in enumerate(names):
            p, q = v[1+2*i], v[2+2*i]
            if p or q: out[nm] = (F(p, v0), F(q, v0))
        print("RELATION:", out)

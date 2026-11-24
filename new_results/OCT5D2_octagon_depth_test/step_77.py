import os, sys, time, pickle, re
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
import mpmath as mp
mp.mp.dps = 430
from fractions import Fraction as F
from flint import fmpz_mat

r2 = mp.sqrt(2); X = mp.log(2); Y = mp.log(1+r2)
def q2f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r2
av = pickle.load(open("atom_vals430.pkl","rb"))
Li5v = {k: mp.mpf(v) for k,v in av['Li5v'].items()}
Li4v = {k: mp.mpf(v) for k,v in av['Li4v'].items()}
Li3v = {k: mp.mpf(v) for k,v in av['Li3v'].items()}
Li2v = {k: mp.mpf(v) for k,v in av['Li2v'].items()}
L32_b1 = mp.mpf(av['L32_b1']); L32_h1 = mp.mpf(av['L32_h1'])
z32 = mp.mpf(av['z32']); K = mp.mpf(av['K'])
z5 = mp.zeta(5); z3 = mp.zeta(3); z2 = mp.pi**2/6
W5 = mp.mpf(pickle.load(open("w5_430.pkl","rb"))['W5'])

data = pickle.load(open("symbol_fit2.pkl","rb"))
cols = data['cols']
xsol = pickle.load(open("symbol_fit2_solution.pkl","rb"))['x']
nz = [(label, kind, xv) for (label, kind), xv in zip(cols, xsol) if xv.p != 0 or xv.q != 0]
def term_value(label):
    m = re.match(r"Li5\((.*)\)$", label)
    if m: return Li5v[m.group(1)]
    m = re.match(r"Li4\((.*)\)\*([XY])$", label)
    if m: return Li4v[m.group(1)] * (X if m.group(2)=='X' else Y)
    for pat, fac in [(r"Li3\((.*)\)\*X\^2/2$", X**2/2), (r"Li3\((.*)\)\*X\*Y$", X*Y), (r"Li3\((.*)\)\*Y\^2/2$", Y**2/2)]:
        m = re.match(pat, label)
        if m: return Li3v[m.group(1)]*fac
    for pat, fac in [(r"Li2\((.*)\)\*X\^3/6$", X**3/6), (r"Li2\((.*)\)\*X\^2Y/2$", X**2*Y/2),
                     (r"Li2\((.*)\)\*XY\^2/2$", X*Y**2/2), (r"Li2\((.*)\)\*Y\^3/6$", Y**3/6)]:
        m = re.match(pat, label)
        if m: return Li2v[m.group(1)]*fac
    m = re.match(r"Li3\((.*)\)\*Li2\((.*)\)$", label)
    if m: return Li3v[m.group(1)] * Li2v[m.group(2)]
    m = re.match(r"Li2\((.*)\)\*Li2\((.*)\)\*([XY])$", label)
    if m: return Li2v[m.group(1)] * Li2v[m.group(2)] * (X if m.group(3)=='X' else Y)
    if label == "Li32(b,1)": return L32_b1
    if label == "Li32(1/sqrt2,1)": return L32_h1
    raise ValueError(label)
cand = sum(q2f(xv)*term_value(label) for label, kind, xv in nz)
R = W5 - cand
print("R =", mp.nstr(R, 25))
invis = {'z5': z5, 'z2*z3': z2*z3, 'z3*X^2': z3*X**2, 'z3*X*Y': z3*X*Y, 'z3*Y^2': z3*Y**2,
         'z2*X^3': z2*X**3, 'z2*X^2Y': z2*X**2*Y, 'z2*XY^2': z2*X*Y**2, 'z2*Y^3': z2*Y**3,
         'z2^2*X': z2**2*X, 'z2^2*Y': z2**2*Y, 'zeta(3,2)': z32, 'K': K}
names = list(invis)
vals = [R] + [invis[n] for n in names]
vv = []
for v in vals: vv += [v, r2*v]
n = len(vv); C = mp.mpf(10)**390
rows = []
for i in range(n):
    row = [0]*(n+1); row[i] = 1; row[n] = int(mp.nint(C*vv[i])); rows.append(row)
M = fmpz_mat(rows).lll()
cands = []
for i in range(n):
    v = [int(M[i,j]) for j in range(n)]
    if v[0] == 0: continue
    res = abs(sum(mp.mpf(v[j])*vv[j] for j in range(n)))
    cands.append((max(abs(xx) for xx in v), res, v))
cands.sort(key=lambda z: z[0])
print("best candidates at dps430 (threshold ~1e15):")
for maxc, res, v in cands[:4]:
    print("  maxc=%.2e res=%s ratio=%s" % (maxc, mp.nstr(res,3), mp.nstr(res*maxc,3)))
    if res < mp.mpf('1e-200'):
        v0 = v[0]; out = {}
        for i, nm in enumerate(names):
            p, q = v[1+2*i], v[2+2*i]
            if p or q: out[nm] = (F(p, v0), F(q, v0))
        print("  RELATION:", out)

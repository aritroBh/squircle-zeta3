import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
from fractions import Fraction as F
import mpmath as mp
from oct5d3_basis import atom_specs, C_SPECS
from oct5d3_atoms import INVISIBLE5, prod, k32_value

sf = pickle.load(open("symfit_m4f.pkl","rb"))
x0, labels = sf['x0'], sf['vis_labels']
Kred = pickle.load(open("kernel_red.pkl","rb"))['Kred']
mixed = pickle.load(open("mixed61_vals.pkl","rb"))
names = pickle.load(open("basis_keep.pkl","rb"))['names']
FORMS = pickle.load(open("closed_forms_lower.pkl","rb"))
rel1500 = {k: tuple(F(x) for x in v) for k, v in pickle.load(open("rel1500_frac.pkl","rb")).items()}
spec_map = {lab: spec for lab, kind, spec in atom_specs()}
vis = [(lab, spec_map[lab]) for lab in labels]
def q2f(z):
    return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*mp.sqrt(2)
def spec_of(label):
    parts = label.split('*'); out = []
    for p in parts:
        if p.startswith('lam^'): out.append(('lam', int(p[4:])))
        elif p == 'lam': out.append(('lam', 1))
        elif p.startswith('L^'): out.append(('L', int(p[2:])))
        elif p == 'L': out.append(('L', 1))
        elif p == 'pi2': out.append(('pi2', 1))
        elif p == 'pi2^2': out.append(('pi2', 2))
        elif p == 'z3': out.append(('z3', 1))
        elif p.startswith('Li2(b)^'): out.append(('Li2(b)', int(p[7:])))
        else: out.append((p, 1))
    return out

D = 800
mp.mp.dps = D
A = {k: mp.mpf(v) for k, v in pickle.load(open("atoms_1650.pkl","rb")).items()}
r2 = mp.sqrt(2)
def cf2(p, q):
    return mp.mpf(p.numerator)/p.denominator + mp.mpf(q.numerator)/q.denominator*r2
# lower weights
tot = mp.mpf(0)
for lvl in (1, 2, 3, 4):
    for label, (a, b) in FORMS[f'Wt{lvl}'].items():
        tot += cf2(F(a), F(b))*prod(A, spec_of(label))
# Wt5: classical + C + REL1500
V = [prod(A, spec) for lab, spec in vis]
classical = sum(q2f(x0[j])*V[j] for j in range(len(V)))
Cval = sum(q2f(cf)*prod(A, spec) for lab, spec, cf in C_SPECS)
kvals = []
for row in Kred:
    v = mp.mpf(0)
    for j in range(len(labels)):
        p, q = row[2*j], row[2*j+1]
        if p or q: v += (mp.mpf(p) + mp.mpf(q)*r2)*V[j]
    kvals.append(v)
allvals = dict(list(zip([nm for nm, sp in INVISIBLE5], [prod(A, sp) for nm, sp in INVISIBLE5]))
               + [("K32", k32_value())]
               + list(zip([f"MIX{i+1}" for i in range(7)], [mp.mpf(val) for val, desc in mixed]))
               + list(zip([f"K{i+1}" for i in range(22)], kvals)))
Wt5val = classical + Cval + sum(cf2(p, q)*allvals[nm] for nm, (p, q) in rel1500.items())
total_over_u = tot + Wt5val
u = mp.sqrt(2+r2)
P5 = u*total_over_u
p5 = pickle.load(open("<project>/research/08_theory/OCT5D_octagon_order5_depth/p5_800.pkl","rb"))
print("closed-form P5/u =", mp.nstr(total_over_u, 40))
print("input P5/u (800) =", mp.nstr(mp.mpf(p5['P5over_u']), 40))
print("difference (P5/u):", mp.nstr(abs(total_over_u - mp.mpf(p5['P5over_u'])), 3))
print("closed-form P5   =", mp.nstr(P5, 40))
print("input P5 (800)   =", mp.nstr(mp.mpf(p5['P5']), 40))
print("difference (P5)  :", mp.nstr(abs(P5 - mp.mpf(p5['P5'])), 3))

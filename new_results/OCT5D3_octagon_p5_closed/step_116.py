import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from fractions import Fraction as F
import mpmath as mp
from oct5d3_basis import atom_specs
from oct5d3_atoms import INVISIBLE5, prod, li32_xy
from qsqrt2 import QS2 as Q

wp = pickle.load(open("wt5_primitive.pkl","rb"))
beta, li32c, planner = wp['beta'], wp['li32'], wp['planner']
spec_map = {lab: spec for lab, kind, spec in atom_specs()}
D = 800
mp.mp.dps = D
A = {k: mp.mpf(v) for k, v in pickle.load(open("atoms_1650.pkl","rb")).items()}
r2 = mp.sqrt(2); b = r2-1
def q2f(z):
    return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r2
LI32V = {
 'Li32(b,1)': li32_xy(b, 1), 'Li32(1/sqrt2,1)': li32_xy(1/r2, 1),
 'Li32(-b,1)': li32_xy(-b, 1), 'Li32(b,-1)': li32_xy(b, -1),
 'Li32(-b,-1)': li32_xy(-b, -1), 'Li32(1/sqrt2,-1)': li32_xy(1/r2, -1),
 'Li32(1-1/sqrt2,1)': li32_xy(1-1/r2, 1), 'Li32(2-sqrt2,1)': li32_xy(2-r2, 1),
 'Li32(-b^2,1)': li32_xy(-b*b, 1),
}
tot = mp.mpf(0)
for lab, cf in beta.items():
    tot += q2f(cf)*prod(A, spec_map[lab])
for lab, cf in li32c.items():
    tot += q2f(cf)*LI32V[lab]
for nm, cf in planner.items():
    spec = dict(INVISIBLE5)[nm] if isinstance(dict(INVISIBLE5), dict) else None
for nm, sp in INVISIBLE5:
    if nm in planner:
        tot += q2f(planner[nm])*prod(A, sp)
Wt5 = mp.mpf(pickle.load(open("wparts_true_1650.pkl","rb"))['Wt5'])
print("primitive expansion value:", mp.nstr(tot, 40))
print("Wt5                    :", mp.nstr(Wt5, 40))
print("difference:", mp.nstr(abs(tot-Wt5), 3))

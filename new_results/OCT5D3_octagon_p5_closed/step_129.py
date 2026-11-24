import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
from oct5d3_basis import atom_list, li32_symbol, add, scale, ZERO, Q
from fractions import Fraction as F

wp = pickle.load(open("wt5_primitive.pkl","rb"))
beta, li32c = wp['beta'], wp['li32']
atoms = atom_list()
sym_map = {a[0]: a[2] for a in atoms if a[1] == 'vis'}
# classical symbol
S_closed = {}
for lab, cf in beta.items():
    if cf != ZERO:
        S_closed = add(S_closed, scale(sym_map[lab], cf))
# Li32 symbol
LI32_ARGS = {'Li32(b,1)': ('b',1), 'Li32(1/sqrt2,1)': ('1/sqrt2',1), 'Li32(-b,1)': ('-b',1),
             'Li32(b,-1)': ('b',-1), 'Li32(-b,-1)': ('-b',-1), 'Li32(1/sqrt2,-1)': ('1/sqrt2',-1),
             'Li32(1-1/sqrt2,1)': ('1-1/sqrt2',1), 'Li32(2-sqrt2,1)': ('2-sqrt2',1), 'Li32(-b^2,1)': ('-b^2',1)}
from oct5d3_basis import B, H, R2, ONE
argmap = {'b': B, '-b': -B, '1/sqrt2': H, '-b^2': -B*B, '1-1/sqrt2': ONE-H, '2-sqrt2': 2-R2}
for lab, cf in li32c.items():
    if cf != ZERO:
        nm, y = LI32_ARGS[lab]
        S_closed = add(S_closed, scale(li32_symbol(argmap[nm], Q(y)), cf))
# target symbol
base2 = "<project>/research/08_theory/OCT5D2_octagon_depth_test/"
raw = pickle.load(open(base2+"block_symbols.pkl","rb"))
W5sym = {tuple(map(int,w)): Q(F(p), F(q)) for w,(p,q) in raw['total'].items()}
diff = add(S_closed, scale(W5sym, -1))
print("S(closed form of Wt5) - S(W5oct) = 0 exactly:", diff == {})
if diff:
    print("NONZERO difference terms:", {k: v for k, v in list(diff.items())[:5]})
from symengine import obstruction
print("obstruction of closed form:", obstruction(S_closed))

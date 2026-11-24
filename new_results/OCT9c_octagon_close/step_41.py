import sys
sys.path.insert(0,'.')
from oct9c_common import *
D = load_all()
engs = {e: get_engine(e) for e in range(4)}
for e in range(4):
    B = D['brackets'][f'B{e}']
    atoms = set()
    for t in sp.Add.make_args(sp.expand(B)):
        for f in sp.Mul.make_args(t):
            b_,e_ = f.as_base_exp()
            atoms.add(b_)
    Ms = sorted([a for a in atoms if str(a).startswith('M_')], key=str)
    print(f"edge {e}: M-atoms:", [str(m) for m in Ms])
    for m in Ms:
        print('     ', m, engs[e].ATOMS[m])
# engine Lam atoms:
for e in range(4):
    G = engs[e]
    print(e, 'Lam1..4:', [G.ATOMS[s] for s in [G.Lam1,G.Lam2,G.Lam3,G.Lam4]])

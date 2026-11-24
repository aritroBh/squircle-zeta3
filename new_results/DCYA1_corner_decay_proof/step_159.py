import sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
from arb_ei import E_I_sc, E_I_cell, iv, arb, alpha
for cell in [(100,121),(100,101),(1600,1626),(1600,1604)]:
    Ncell=iv(str(cell[0]),str(cell[1]))
    try:
        v=E_I_sc(Ncell)
        print(f"{cell}: {float(v):.3f}")
    except Exception as e:
        print(f"{cell}: ERR {type(e).__name__} {e}")

import sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
from arb_ei import E_I_sc, iv
for cell in [(1410,1415),(1410,1420),(1318,1325),(1310,1318),(1384,1390),(1600,1604)]:
    try:
        v=float(E_I_sc(iv(str(cell[0]),str(cell[1]))))
        print(f"{cell}: E_I/sc={v:.3f}")
    except Exception as e:
        print(f"{cell}: ERR {type(e).__name__}")

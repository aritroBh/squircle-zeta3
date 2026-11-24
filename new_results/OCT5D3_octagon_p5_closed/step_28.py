import os
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
md = """# OCT5D3 — closed form of the regular octagon's P⁽⁵⁾(1) — STATUS: IN PROGRESS

Inputs: OCT5D2 (depth test, engines, symbol fit), OCT5D (parsed_mons, p5_800),
M4e/M4f (squircle W5 closed form: basis + conventions).

## Plan
1. Classical representative of W5oct/u in the M4f basis (14 Li5 args + M4f products + C), exact symbol algebra.
2. G-table to >=1300 digits; W5 and lower-weight parts; cross-check p5_800.
3. Closure of R against the completed invisible basis (planner diagnosis), two precisions.
4. Full assembly of P⁽⁵⁾(1); zeta(5), zeta(2)zeta(3) coefficients vs prop:order5 (m=4).

(progress log below)
"""
open("result.md","w").write(md)
print("wrote result.md")
# Check block_symbols.pkl format
import sys, pickle
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
sys.path.insert(0, "<project>/research/08_theory/OCT5D2_octagon_depth_test")
base2 = "<project>/research/08_theory/OCT5D2_octagon_depth_test/"
bs = pickle.load(open(base2+"block_symbols.pkl","rb"))
print(bs.keys())
print({k: str(v)[:150] for k,v in bs.items()})

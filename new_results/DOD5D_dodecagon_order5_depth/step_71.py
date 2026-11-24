import sys, json, pickle
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
from fractions import Fraction as F
from qsqrt2 import QS2 as Q
import obstruction32 as OB

# ---- M4c control ----
symj = json.loads(open("<project>/research/08_theory/M4c_w5_closed_form/symbol.json").read())['symbol']
Sm4 = {tuple(map(int, w)): Q(F(p), F(q)) for w, (p, q) in symj.items()}
table, (n3, n2) = OB.obstruction(Sm4, 2)
print("M4c W5, Lie3 names:", n3, "Lie2:", n2)
for i, row in enumerate(table):
    print("  ", n3[i], [str(v) for v in row])
print("expect [X,[X,Y]]x[X,Y] = 84, [Y,[X,Y]]x[X,Y] = -38*sqrt2:", Q(84), Q(0,-38))

# ---- octagon control ----
bs = pickle.load(open("<project>/research/08_theory/OCT5D2_octagon_depth_test/block_symbols.pkl","rb"))
def load(d): return {tuple(map(int, w)): Q(F(p), F(q)) for w, (p, q) in d.items()}
for blk in ['local', 'single', 'mixed', 'total']:
    S = load(bs['blocks'][blk] if blk != 'total' else bs['total'])
    table, _ = OB.obstruction(S, 2)
    print("octagon", blk, ":", [str(row[0]) for row in table])
print("expect total (-1485+1035 sqrt2, 900-540 sqrt2):", Q(-1485,1035), Q(900,-540))

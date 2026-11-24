
import os, json
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
os.chdir(WD)
# E_{5,m} coefficient tables
e = json.load(open('e5m_binomial_coeffs.json'))
KEYS5, POLES5 = ["50","41","32"], ["0","a","-a","-b"]
from fractions import Fraction as F
def fmt(rat, sq):
    r, s = F(rat), F(sq)
    if s == 0: return str(r)
    if r == 0: return f"{s}*sqrt2"
    return f"({r}{'+' if s>0 else '-'}{abs(s)}*sqrt2)"
for m in range(1,6):
    terms = []
    for key in KEYS5:
        i,j = int(key[0]), int(key[1])
        for pole in POLES5:
            cell = e[f"{key}|{pole}"]
            if F(cell['rat'][m-1]) or F(cell['sqrt2'][m-1]):
                terms.append(f"{fmt(cell['rat'][m-1], cell['sqrt2'][m-1])} N({i},{j},{pole})")
    print(f"E5,{m}: " + " + ".join(terms))
    print()

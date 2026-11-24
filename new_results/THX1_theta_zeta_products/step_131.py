
import pickle
from fractions import Fraction as Fr
reds9 = pickle.load(open('I_reductions_w9.pkl','rb'))
basis_names = pickle.load(open('I_reductions_w89.pkl','rb'))['basis9']
def fmt(c):
    if c == 1: return ""
    return f"{c}"
for (m, j), coeffs in sorted(reds9.items()):
    terms = []
    for c, nm in zip(coeffs, basis_names):
        if c:
            terms.append(f"{c}*{nm}" if c != 1 else nm)
    print(f"I_{{{m},{j}}} = " + " + ".join(terms))

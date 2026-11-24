import pickle, sys
sys.path.insert(0, ".")
from fractions import Fraction as F
from qs3_symbol import Q3V
vecdata = pickle.load(open("mpl_cobrackets.pkl","rb"))
sol = pickle.load(open("li32_solution_raw.pkl","rb"))
# pair up (2i, 2i+1) into QS3 coefficients
terms = {}
for j, xj in sol:
    i = j // 2
    terms.setdefault(i, [F(0), F(0)])
    terms[i][j % 2] = F(xj)
def fmtq3(p, q):
    if q == 0: return str(p)
    if p == 0: return f"{q}r3"
    return f"{p}{q:+}r3"
def fmtz(z):
    p, q = F(z[0]), F(z[1])
    return f"({fmtq3(p,q)})"
out = []
for i, (p, q) in sorted(terms.items()):
    (ab, x, y, _) = vecdata[i]
    out.append((ab, x, y, Q3V(p, q)))
for (ab, x, y, c) in out:
    print(f"  ({fmtq3(c[0],c[1])}) * Li_{ab[0]},{ab[1]}( {fmtz(x)}, {fmtz(y)} )")
print("total terms:", len(out))

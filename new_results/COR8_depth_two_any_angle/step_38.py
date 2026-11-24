
import pickle, numpy as np
from fractions import Fraction as Fr
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve_data.pkl','rb'))
R, li4, prods, comps = D['R'], D['li4'], D['prods'], D['comps']
# word index per exponent: symbols are per-exponent e=0 for ansatz; R has e in {-1,0,1}
words = sorted({w for S in [R]+list(li4.values())+list(prods.values()) for w in S}, key=str)
widx = {w:i for i,w in enumerate(words)}
NW = len(words)
print("words:", NW)
# dense vectors mod p
p = (1<<31)-1
def vecs_modp(S):
    v = np.zeros(NW, dtype=np.int64)
    for w, cf in S.items():
        if 0 in cf:
            v[widx[w]] = int(cf[0] % p)
    return v
# incremental basis of products mod p
items = list(prods.items())
basis_idx = []
pivots = []   # word positions
B_rows = []   # basis vectors (mod p), in row-echelon form
for k, (name, S) in enumerate(items):
    v = vecs_modp(S)
    for br, pw in zip(B_rows, pivots):
        if v[pw]:
            v = (v - int(v[pw]) * br) % p
    nz = np.nonzero(v)[0]
    if len(nz):
        pw = int(nz[0])
        inv = pow(int(v[pw]), p-2, p)
        br = (v * inv) % p
        B_rows.append(br); pivots.append(pw); basis_idx.append(k)
print("product space rank (mod p):", len(B_rows))
pivot_words = [words[i] for i in pivots]
basis_names = [items[k][0] for k in basis_idx]
for nm, pw in zip(basis_names, pivot_words):
    print("  ", nm)
pickle.dump({'words':words,'pivot_words':pivot_words,'basis_names':basis_names}, open('prod_basis.pkl','wb'))

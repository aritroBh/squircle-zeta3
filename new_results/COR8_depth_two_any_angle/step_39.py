
import pickle, numpy as np
from fractions import Fraction as Fr
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve_data.pkl','rb'))
R, li4, prods = D['R'], D['li4'], D['prods']
# basis words: tuples of 4 elements of {0,1,2} = {c,m,p}
import itertools
bwords = list(itertools.product(range(3), repeat=4))
bidx = {w:i for i,w in enumerate(bwords)}
def expand_letter(l):
    # l is exponent vector over (c,m,p): returns list of (basis_index, Fr)
    return [(j, l[j]) for j in range(3) if l[j]]
def to81(S, e=0):
    """symbol -> 81-dim Fraction vector at exponent e"""
    v = [Fr(0)]*81
    for w, cf in S.items():
        if e not in cf: continue
        coeff = cf[e]
        # expand multilinearly
        terms = [((), Fr(1))]
        for l in w:
            opts = expand_letter(l)
            terms = [(p+(j,), c0*m) for p, c0 in terms for j, m in opts]
        for p, c0 in terms:
            v[bidx[p]] += coeff*c0
    return v
# convert everything once
li4v = {k: to81(S) for k, S in li4.items()}
prodv = {k: to81(S) for k, S in prods.items()}
Rv = {E: to81(R, E) for E in (-1,0,1)}
# ranks mod p
p = (1<<31)-1
def rank_modp(vecs):
    A = np.array([[x % p for x in v] for v in vecs], dtype=np.int64)
    # simple gaussian elimination
    A = A.copy(); m, n = A.shape; r = 0; cols = []
    for c in range(n):
        piv = -1
        for i in range(r, m):
            if A[i,c] % p: piv = i; break
        if piv < 0: continue
        A[[r,piv]] = A[[piv,r]]
        inv = pow(int(A[r,c]%p), p-2, p)
        A[r] = (A[r]*inv) % p
        for i in range(m):
            if i != r and A[i,c]%p:
                A[i] = (A[i] - A[i,c]*A[r]) % p
        cols.append(c); r += 1
        if r == m: break
    return r, cols
allv = list(li4v.values()) + list(prodv.values())
rtot, _ = rank_modp(allv)
rprod, _ = rank_modp(list(prodv.values()))
rli4, _ = rank_modp(list(li4v.values()))
print("rank Li4:", rli4, " rank products:", rprod, " rank all:", rtot, " (81 max)")
for E in (-1,0,1):
    re, _ = rank_modp(list(li4v.values())+list(prodv.values())+[Rv[E]])
    print(f"rank with R_{E}:", re)
pickle.dump({'li4v':li4v,'prodv':prodv,'Rv':Rv,'bwords':bwords}, open('solve81.pkl','wb'))

add = '''

# ---- weight-5 quotient T^5/D^5 (free Lie_5 coords), NF coefficients ----
_T5 = None
def t5d5():
    global _T5
    if _T5 is None:
        from itertools import product as _ip
        B5 = list(_ip(range(2), repeat=5))
        rows = []
        for lu in range(1, 5):
            for u in _ip(range(2), repeat=lu):
                for v in _ip(range(2), repeat=5 - lu):
                    rows.append(sshuffle({u: Fr(1)}, {v: Fr(1)}))
        Dmat5, piv5 = rref([[s.get(w, Fr(0)) for w in B5] for s in rows])
        free5 = [j for j in range(len(B5)) if j not in piv5]
        _T5 = (B5, Dmat5, piv5, free5)
    return _T5

def quotient5(SE, dimF=2):
    """SE {index-word -> NF} -> tuple of 6 NF coords in T^5/D^5."""
    B5, Dmat5, piv5, free5 = t5d5()
    v = [tuple(SE.get(w, (Fr(0),) * dimF)) for w in B5]
    for i, pj in enumerate(piv5):
        if any(v[pj]):
            q = v[pj]
            v = [tuple(vv[t] - q[t] * Dmat5[i][j] for t in range(dimF)) for j, vv in enumerate(v)]
    for j in range(len(B5)):
        if j not in free5:
            assert not any(v[j]), "not integrable"
    return tuple(v[j] for j in free5)

def li5_sym_nf(P, z):
    """S(Li5(z)) = -[1-z][z]^4, exponent-vector letters, NF coeff."""
    F = P.F
    vz = P.vec(z); v1z = P.vec(F.sub(F.one(), z))
    return {(v1z, vz, vz, vz, vz): (Fr(-1),) + (Fr(0),) * (F.d - 1)}
'''
src = open('gold5_lib.py').read()
open('gold5_lib.py','w').write(src + add)
print("added")
import sys
sys.path.insert(0, '<project>/research/08_theory/GOLD5_golden_corner_w5')
import importlib
import gold5_lib as G
importlib.reload(G)
from fractions import Fraction as Fr
# class of Omega_5 at both corners
for tag in ('gol1', 'gol3'):
    P = G.golden_point(tag)
    S5 = P.omega5()
    SE = G.ft_expand_sym_nf(S5, 2)
    cl = G.quotient5(SE)
    print(f"{tag}: class(Omega_5) = {[G.nf_str(q) for q in cl]}")


# patch engine: linear maps for relabel/specialize (vectors are formal sums)
patch = r'''

# ---- linear-map versions (letters may be formal sums of basis vectors) ----
def _apply_lin(basis_images, v):
    out = [0]*len(next(iter(basis_images.values())))
    basis = [XV,YV,UV,VV,WV] if len(v)==5 else None
    if basis is None:
        # 4-dim target side identity
        basis = [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)]
    for i,b in enumerate(basis):
        img = basis_images.get(b, b)
        for j in range(len(out)):
            out[j] += v[i]*img[j]
    return tuple(out)

def relabel_lin(sym, basis_images):
    out = {}
    for k,c in sym.items():
        nk = tuple(_apply_lin(basis_images, v) for v in k)
        out[nk] = out.get(nk, Fr(0)) + c
    return {k:v for k,v in out.items() if v}

def specialize_lin(sym5, smap):
    return relabel_lin(sym5, smap)
'''
open('cor8b_sym.py','a').write(patch)
import sys; sys.path.insert(0,'.')
import importlib, cor8b_sym
importlib.reload(cor8b_sym)
from cor8b_sym import *

S22 = S_Li22(); S31 = S_Li31(); S13 = S_Li13()
swap = {XV:YV, YV:XV, UV:VV, VV:UV, WV:WV}
def shuffle2(symA, symB):
    out = {}
    for ka, ca in symA.items():
        for kb, cb in symB.items():
            a1,a2 = ka; b1,b2 = kb
            for perm in [(a1,a2,b1,b2),(a1,b1,a2,b2),(a1,b1,b2,a2),(b1,a1,a2,b2),(b1,a1,b2,a2),(b1,b2,a1,a2)]:
                out[perm] = out.get(perm, Fr(0)) + ca*cb
    return {k:v for k,v in out.items() if v}
lhs = sym_add(sym_add(S22, relabel_lin(S22, swap)), S_Li4z('xy'))
rhs = shuffle2(S_Li2('x'), S_Li2('y'))
print("shuffle Li22 residual terms:", len(sym_add(lhs, sym_scale(rhs,-1))))

S3x = { (UV,XV,XV): Fr(-1) }; S1y = { (VV,): Fr(-1) }
def shuffle31(symA, symB):
    out = {}
    for ka, ca in symA.items():
        for kb, cb in symB.items():
            a1,a2,a3 = ka; (b1,) = kb
            for perm in [(a1,a2,a3,b1),(a1,a2,b1,a3),(a1,b1,a2,a3),(b1,a1,a2,a3)]:
                out[perm] = out.get(perm, Fr(0)) + ca*cb
    return {k:v for k,v in out.items() if v}
lhs2 = sym_add(sym_add(S31, relabel_lin(S13, swap)), S_Li4z('xy'))
rhs2 = shuffle31(S3x, S1y)
print("stuffle Li31 residual terms:", len(sym_add(lhs2, sym_scale(rhs2,-1))))

for sgn,name in [(-1,'y=1-c'),(1,'y=1+c')]:
    sm = spec_map(sgn)
    d22 = delta22(specialize_lin(S22, sm))
    d31 = delta22(specialize_lin(S31, sm))
    d13 = delta22(specialize_lin(S13, sm))
    print(f"--- {name}")
    print(" d22(Li22) =", show_d22(d22))
    print(" d22(Li31) =", show_d22(d31))
    print(" d22(Li13) =", show_d22(d13))

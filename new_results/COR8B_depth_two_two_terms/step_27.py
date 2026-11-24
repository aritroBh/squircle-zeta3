
import sys
sys.path.insert(0,'.')
from cor8b_sym import *

# ---------- two-variable symbols ----------
S22 = S_Li22(); S31 = S_Li31(); S13 = S_Li13(); S12 = S_Li12(); S21 = S_Li21(); S11 = S_Li11()
print("S(Li22) terms:", len(S22), " S(Li31):", len(S31), " S(Li13):", len(S13), " S(Li12):", len(S12), " S(Li21):", len(S21), " S(Li11):", len(S11))

# ---------- shuffle identity: Li22(x,y)+Li22(y,x)+Li4(xy) = Li2(x)Li2(y) ----------
# swap x<->y: X<->Y, U<->V, W->W
swap = {XV:YV, YV:XV, UV:VV, VV:UV, WV:WV}
def relabel(sym, mp):
    return {tuple(mp[v] for v in k): c for k,c in sym.items()}
def shuffle2(symA, symB):
    # shuffle product of two weight-2 symbols
    out = {}
    for ka, ca in symA.items():
        for kb, cb in symB.items():
            a1,a2 = ka; b1,b2 = kb
            for perm in [(a1,a2,b1,b2),(a1,b1,a2,b2),(a1,b1,b2,a2),(b1,a1,a2,b2),(b1,a1,b2,a2),(b1,b2,a1,a2)]:
                out[perm] = out.get(perm, Fr(0)) + ca*cb
    return {k:v for k,v in out.items() if v}

lhs = sym_add(sym_add(S22, relabel(S22, swap)), S_Li4z('xy'))
rhs = shuffle2(S_Li2('x'), S_Li2('y'))
diff = sym_add(lhs, sym_scale(rhs, -1))
print("shuffle Li22 residual terms:", len(diff))

# stuffle: Li3(x)Li1(y) = Li31(x,y)+Li13(y,x)+Li4(xy)
# S(Li3(x)) = -(U,X,X); S(Li1(y)) = -(V,)
S3x = { (UV,XV,XV): Fr(-1) }
S1y = { (VV,): Fr(-1) }
def shuffle31(symA, symB):
    out = {}
    for ka, ca in symA.items():
        for kb, cb in symB.items():
            a1,a2,a3 = ka; (b1,) = kb
            for perm in [(a1,a2,a3,b1),(a1,a2,b1,a3),(a1,b1,a2,a3),(b1,a1,a2,a3)]:
                out[perm] = out.get(perm, Fr(0)) + ca*cb
    return {k:v for k,v in out.items() if v}
lhs2 = sym_add(sym_add(S31, relabel(S13, swap)), S_Li4z('xy'))
rhs2 = shuffle31(S3x, S1y)
diff2 = sym_add(lhs2, sym_scale(rhs2, -1))
print("stuffle Li31 residual terms:", len(diff2))

# ---------- delta22 at two-variable level should be checked on specialization ----------
# specialize to y=1-c, y=1+c and compute delta22
for sgn,name in [(-1,'y=1-c'),(1,'y=1+c')]:
    sm = spec_map(sgn)
    d22 = delta22(specialize(S22, sm))
    d31 = delta22(specialize(S31, sm))
    d13 = delta22(specialize(S13, sm))
    print(f"--- {name}")
    print(" d22(Li22) =", show_d22(d22))
    print(" d22(Li31) =", show_d22(d31))
    print(" d22(Li13) =", show_d22(d13))

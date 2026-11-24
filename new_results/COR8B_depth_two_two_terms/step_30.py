
import sys; sys.path.insert(0,'.')
from cor8b_sym import *
S22 = S_Li22()
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
diff = sym_add(lhs, sym_scale(rhs,-1))
B5 = ['X','Y','U','V','W']
def show5(sym):
    for k,c in sorted(sym.items()):
        print("   ", c, "*", "["+"][".join(B5[list(v).index(1)] if v.count(1)==1 and v.count(0)==len(v)-1 else str(v) for v in k)+"]")
print("RESIDUAL:")
show5(diff)
print("\nS(Li22) two-variable:")
show5(S22)

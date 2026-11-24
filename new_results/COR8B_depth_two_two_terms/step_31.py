
import sys; sys.path.insert(0,'.')
from cor8b_sym import *
import itertools
B5 = [XV,YV,UV,VV,WV]
def normalize5(sym):
    # expand every tensor multilinearly into basis monomials
    out = {}
    for k,c in sym.items():
        # for each slot, get list of (basis_index, multiplicity)
        opts = []
        for v in k:
            o = [(i, vi) for i,vi in enumerate(v) if vi]
            opts.append(o)
        for combo in itertools.product(*opts):
            mult = Fr(c)
            for i,vi in combo: mult *= vi
            nk = tuple(B5[i] for i,vi in combo)
            out[nk] = out.get(nk, Fr(0)) + mult
            if out[nk]==0: del out[nk]
    return out

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
diff = normalize5(sym_add(lhs, sym_scale(rhs,-1)))
print("shuffle Li22 (normalized) residual terms:", len(diff), diff if len(diff)<5 else '')

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
diff2 = normalize5(sym_add(lhs2, sym_scale(rhs2,-1)))
print("stuffle Li31 (normalized) residual terms:", len(diff2), diff2 if len(diff2)<5 else '')

# also: Li22(x,y)+Li22(y,x)+Li4(xy)=Li2(x)Li2(y) implies after specialization numeric check already done (step_01 residual 1e-53)
# save normalized two-variable symbols for the record
import json
def sym_to_json(sym):
    return {"|".join("".join(map(str,v)) for v in k): str(c) for k,c in sym.items()}
json.dump({'Li22':sym_to_json(S_Li22()),'Li31':sym_to_json(S_Li31()),'Li13':sym_to_json(S_Li13()),
           'Li12':sym_to_json(S_Li12()),'Li21':sym_to_json(S_Li21()),'Li11':sym_to_json(S_Li11())},
          open('twovar_symbols.json','w'), indent=1)
print("saved twovar_symbols.json")

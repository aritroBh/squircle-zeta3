
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ
import math

e_rm = pickle.load(open('e_rm_ext.pkl','rb'))
def E(r, m): return Fr(e_rm.get(f"{r},{m}", '0'))
def I_word_combo(m, j):
    return {(0,)*(m+1) + (-1,)*j: Fr(math.factorial(j)*math.factorial(m))}
def Dtilde_combo(n, m):
    out = {}
    for j in range(1, n-2*m):
        e = E(m, n-1-j)
        if e == 0: continue
        c = 2*e/math.factorial(j)
        for w, vv in I_word_combo(n-1-j, j).items():
            out[w] = out.get(w, Fr(0)) + c*vv
    return {k: v for k, v in out.items() if v}
def zeta_word(n): return {(0,)*(n-1)+(1,): Fr(-1)}
def zprod(*ns):
    out = {(): Fr(1)}
    for n in ns: out = rul1.sh_mul(out, zeta_word(n))
    return {k: v for k, v in out.items() if v}
Dtarget = pickle.load(open('Dtilde_relations.pkl','rb'))
def V_target(n, m):
    out = Dtilde_combo(n, m)
    coeffs = [Fr(c) for c in Dtarget[f"{n},{m}"]]
    for c, zb in zip(coeffs, [zprod(7), zprod(2,5), zprod(2,2,3)]):
        for w, vv in zb.items():
            out[w] = out.get(w, Fr(0)) - c*vv
    return {k: v for k, v in out.items() if v}

rels_data = pickle.load(open('rels7_dsdual.pkl','rb'))
rels = [{k: Fr(n, d) for k, (n, d) in r.items()} for fam, tag, r in rels_data]
cols7 = sorted(rul1.conv_words_of_weight(7))
st = pickle.load(open('struct_w7.pkl','rb'))
B = st['indep_rows']; colpos = st['colpos']; pivcols = st['pivcols']
r = len(B)
posmap = {c: jj for jj, c in enumerate(pivcols)}

t=time.time()
Sq = {}
for i, bi in enumerate(B):
    for k, v in rels[bi].items():
        j = posmap.get(colpos[k])
        if j is not None:
            Sq.setdefault(j, {})[i] = QQ(v.numerator, v.denominator)
MS = DomainMatrix(Sq, (r, r), QQ).to_sparse()
V = V_target(7, 0)
bdata = {}
for k, v in V.items():
    j = posmap.get(colpos[k])
    if j is not None:
        bdata[j] = {0: QQ(v.numerator, v.denominator)}
b = DomainMatrix(bdata, (r, 1), QQ).to_sparse()
print("built system:", time.time()-t)
t=time.time()
sol = MS.transpose().lu_solve(b)
print("lu_solve time:", time.time()-t)
solD = sol.to_dict_sympy()
nz = sum(1 for i in range(r) if solD.get((i,0),0) != 0)
print("nonzero solution entries:", nz)
# verify exactly: h . R_B == V
acc = {}
for i in range(r):
    hi = solD.get((i,0),0)
    if hi == 0: continue
    hF = Fr(hi.p, hi.q)
    for k, v in rels[B[i]].items():
        acc[k] = acc.get(k, Fr(0)) + hF*v
for k, v in V.items():
    acc[k] = acc.get(k, Fr(0)) - v
acc = {k: v for k, v in acc.items() if v}
print("EXACT VERIFICATION:", "PASS" if acc == {} else f"FAIL ({len(acc)} bad cols)")
if acc == {}:
    pickle.dump({i: str(solD[(i,0)]) for i in range(r) if solD.get((i,0),0)!=0}, open('cert_D7_0.pkl','wb'))
    open('progress.log','a').write(time.strftime('%H:%M:%S ') + "D~_{7,0} PROVED via exact square solve + verification\n")

import cor8r_sym as S
import sympy as sp
from fractions import Fraction
import itertools
cc = S.cc
A = 1/(1+cc); B = 2/(1+cc)
L = 3
values = [A, B]
unk = S._bad_words(L, values)
idx = {w:i for i,w in enumerate(unk)}
alphabet = [sp.Integer(0)]+values
rows = []
for lu in range(1,L):
    lv = L-lu
    for u in itertools.product(alphabet, repeat=lu):
        su = S.symG(u)
        for v in itertools.product(alphabet, repeat=lv):
            sv = S.symG(v)
            prod = S.S_shuffle(su,sv)
            row = [Fraction(0)]*len(unk)
            rhs = dict(prod)
            ok = False
            for w in S._shuffle_words(u,v):
                if w in idx:
                    row[idx[w]] += 1; ok = True
                else:
                    sw = S.symG(w)
                    rhs = S.S_add(rhs, {ww:S.LP_neg(cf) for ww,cf in sw.items()})
            if ok:
                rows.append((u,v,row,rhs))
# solve with the used-rows elimination (no consistency raise)
M = [(list(r[2]), r[3], (r[0],r[1])) for r in rows]
sol = {}
R = 0
for col in range(len(unk)):
    piv = None
    for r in range(R, len(M)):
        if M[r][0][col] != 0: piv = r; break
    if piv is None: continue
    M[R], M[piv] = M[piv], M[R]
    inv = Fraction(1)/M[R][0][col]
    M[R] = ([x*inv for x in M[R][0]], S.S_scaleLP_frac(M[R][1], inv), M[R][2])
    for r in range(len(M)):
        if r!=R and M[r][0][col]!=0:
            f = M[r][0][col]
            M[r] = ([a-f*b for a,b in zip(M[r][0],M[R][0])], S.S_add(M[r][1], S.S_scaleLP_frac(M[R][1],-f)), M[r][2])
    R += 1
    if R==len(M): break
# solution: pivot rows
for rr in range(R):
    row, s, uv = M[rr]
    pc = [i for i,x in enumerate(row) if x!=0]
    if len(pc)==1:
        sol[pc[0]] = s
print("solved unknowns:", len(sol))
# residuals of original rows
def resid(row, rhs):
    r = dict(rhs)
    for i,c in enumerate(row):
        if c and i in sol:
            r = S.S_add(r, S.S_scaleLP_frac(sol[i], -c))
        elif c and i not in sol:
            return "UNSOLVED"
    return r
bad = 0
for u,v,row,rhs in rows:
    r = resid(row,rhs)
    if r and r!="UNSOLVED":
        bad += 1
        if bad<=3:
            print("\nRESID u,v=",u,v)
            for w,cf in r.items():
                print("   ", w, cf)
print("bad rows:", bad)

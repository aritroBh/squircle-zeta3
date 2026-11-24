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
# solve via sympy RREF with symbolic rhs: solve for each unknown's symbol using first indep rows
# Build augmented matrix where rhs is a fresh symbol per row? Instead: solve numerically on random c specialization.
# Simpler: find 16 independent rows, solve, then test all.
import copy
M = [(list(r[2]), r[3], (r[0],r[1])) for r in rows]
R = 0; used = []
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
    used.append(M[R][2]); R += 1
    if R==len(M): break
print("used rows:", used)
print("\nFailing original equations (0-row with nonzero rhs):")
for row,s,uv in M:
    if all(x==0 for x in row) and s:
        print("  u,v =", uv, " rhs #terms:", len(s))

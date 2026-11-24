import cor8r_sym as S
import sympy as sp
from fractions import Fraction
import itertools
cc = S.cc
TWO = sp.Integer(2); B = 2/(1+cc)
for L, values in [(2,[TWO,B]),(3,[TWO,B]),(4,[TWO,B])]:
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
    M = [(list(r[2]), r[3], (r[0],r[1])) for r in rows]
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
    sol = {}
    for rr in range(R):
        row, s, uv = M[rr]
        pc = [i for i,x in enumerate(row) if x!=0]
        if len(pc)==1: sol[pc[0]] = s
    bad = 0
    for u,v,row,rhs in rows:
        r = dict(rhs)
        for i,cf in enumerate(row):
            if cf: r = S.S_add(r, S.S_scaleLP_frac(sol[i], -cf))
        if r: bad += 1
    print(f"L={L}: unknowns {len(unk)}, rows {len(rows)}, solved {len(sol)}, inconsistent rows {bad}")

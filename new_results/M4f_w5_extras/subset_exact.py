"""M4f step 3: exact-coordinate subset search.

Value space V = Q^86 / LAM, LAM = 14 exact relation vectors (from overcomplete_1200.json,
verified to ~1300 digits). Basis B = 72 kept columns. Columns 0..34 kv, 35..69 s2kv,
70 L, 71 s2L, 72..78 u0..u6, 79..85 s2u0..s2u6.
Atom slots 0..64 p / 65..129 q (65 S0 atoms), 130/131 Lchi p/q, 132..138 Li5(z_i) p, 139..145 q.
For every subset S of the 7 extras: exact solve M_S x = tbar; if closing, optimize
t in x0 + N_S to minimise (max|num|, lcm den, terms) of the atom expansion
(after Li5(b^2)->16(Li5(b)+Li5(-b)), Li5(-1)->-15/16 zeta(5)).
"""
import json, sys, itertools, math
from fractions import Fraction as F
from math import gcd, lcm
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
import nice_search as NS

def log(s): print(s, flush=True)

def rref_solve(rows, ncols, nrhs):
    mat, rank = flint.fmpq_mat(rows).rref()
    pivots = [next((j for j in range(mat.ncols()) if mat[i,j]), None) for i in range(rank)]
    return mat, rank, pivots

def main():
    oc = json.loads((HERE/'overcomplete_1200.json').read_text())
    deps86 = []
    for dep in oc['deps']:
        assert len(dep['kept']) == 86 and dep['vec'][0] == 0
        deps86.append([F(c) for c in dep['vec'][1:]])
    assert len(deps86) == 14
    rel = oc['target_rels'][0]
    c0 = F(rel['vec'][0]); relv = [F(c) for c in rel['vec'][1:]]
    assert len(relv) == 72
    kept = oc['kept']
    kidx = []
    for nm in kept:
        if nm == 'L': kidx.append(70)
        elif nm == 's2L': kidx.append(71)
        elif nm.startswith('s2kv'): kidx.append(35+int(nm[4:]))
        elif nm.startswith('kv'): kidx.append(int(nm[2:]))
        elif nm.startswith('s2u'): kidx.append(79+int(nm[3:]))
        elif nm.startswith('u'): kidx.append(72+int(nm[1:]))
        else: raise ValueError(nm)
    # coordinate solve: matrix [B | D] (86 x 86), rhs = e_c for all 86 columns c
    Bcols = [[F(int(i==j)) for i in range(86)] for j in kidx]   # columns of B as 86-vectors
    M = [[Bcols[j][i] for j in range(72)] + [deps86[d][i] for d in range(14)] for i in range(86)]
    rhs = [[F(int(i==j)) for j in range(86)] for i in range(86)]
    rows = [[str(c) for c in M[i] + rhs[i]] for i in range(86)]
    mat, rank = flint.fmpq_mat(rows).rref()
    assert rank == 86, f'[B|D] not invertible: rank {rank}'
    inv = [[F(str(mat[i,86+j])) for j in range(86)] for i in range(86)]  # inverse of [B|D]
    # coords in B of column c = first 72 entries of inv * e_c = column c of inv, first 72 rows
    coords = [[inv[i][c] for i in range(72)] for c in range(86)]
    # check: dropped u columns reconstruct as combos; print readable relations
    names_u = ['u0(=-b/2)','u1(=2b)','u2((2-s2)/4)','u3((2+s2)/4)','u4(b^4)','u5((4-3s2)/8)','u6(1-b^4)']
    tbar = [F(-relv[i], c0) for i in range(72)]
    # exact u-relations in terms of basis columns
    for c in list(range(72,79)):
        if c-72 in [k-72 for k in kidx if k>=72]: continue
        cv = coords[c]
        sup = [(kidx[i], cv[i]) for i in range(72) if cv[i]]
        usup = [(names_u[j-72] if j>=72 else ('L' if j==70 else 's2L' if j==71 else f'kv{j}' if j<35 else f's2kv{j-35}'), v)
                for j,v in sup]
        log(f'exact relation: {names_u[c-72]} = {usup}')
    # atom map
    sol = NS.solve(list(NS.S0), 'subset')
    names = sol['names']; xp0, xq0 = sol['xp'], sol['xq']
    K = json.loads((HERE/'nice_ksat_s0.json').read_text())['K']
    pool = json.loads((HERE/'pool.json').read_text())['pool']
    s0 = set(NS.S0)
    extras = [e for e in pool if NS.W.Q(F(e['p']),F(e['q'])) not in s0]
    assert len(extras) == 7
    NATOM = 146
    def acol(c):
        v = [F(0)]*NATOM
        if c < 35:
            for i,a in enumerate(K[c]): v[i] = F(a)
        elif c < 70:
            for i,a in enumerate(K[c-35]): v[65+i] = F(a)
        elif c == 70: v[130] = F(1)
        elif c == 71: v[131] = F(1)
        elif c < 79:
            i = c-72; e = extras[i]
            v[132+i] = F(1)
            for j,a in enumerate(e['xp']): v[j] -= F(a)
            for j,a in enumerate(e['xq']): v[65+j] -= F(a)
        else:
            i = c-79; e = extras[i]
            v[139+i] = F(1)
            for j,a in enumerate(e['xq']): v[j] -= 2*F(a)
            for j,a in enumerate(e['xp']): v[65+j] -= F(a)
        return v
    ACOLS = [acol(c) for c in range(86)]
    vpart = [F(0)]*NATOM
    for i,a in enumerate(xp0): vpart[i] = a
    for i,a in enumerate(xq0): vpart[65+i] = a
    # reduction indices
    ib2 = names.index('Li5((3-2*sqrt2))'); ibp = names.index('Li5((-1+1*sqrt2))')
    ibm = names.index('Li5((1-1*sqrt2))'); im1 = names.index('Li5(-1)'); iz5 = names.index('z5')
    def reduced(v):
        v = list(v)
        for s, off in ((0,0),(65,65)):
            a = v[ib2+off]
            if a:
                v[ibp+off] += 16*a; v[ibm+off] += 16*a; v[ib2+off] = F(0)
            a = v[im1+off]
            if a:
                v[iz5+off] += F(-15,16)*a; v[im1+off] = F(0)
        return v
    def height(v):
        v = reduced(v)
        nz = [a for a in v if a]
        mx = max(abs(a.numerator) for a in nz)
        dn = lcm(*[a.denominator for a in nz])
        return mx, dn, len(nz)
    # subsets
    results = []
    uall = list(range(7))
    for size in range(3, 8):
        for S in itertools.combinations(uall, size):
            allow = list(range(72)) + [72+i for i in S] + [79+i for i in S]
            ncol = len(allow)
            Mx = [[str(coords[c][i]) for c in allow] + [str(tbar[i])] for i in range(72)]
            mat, rank = flint.fmpq_mat(Mx).rref()
            pivots = [next((j for j in range(mat.ncols()) if mat[i,j]), None) for i in range(rank)]
            if any(p == ncol for p in pivots):   # tbar not in span
                continue
            free = [j for j in range(ncol) if j not in pivots]
            x0 = [F(0)]*ncol
            for i,p in enumerate(pivots):
                if p is not None and p < ncol: x0[p] = F(str(mat[i,ncol]))
            ns = []
            for f in free:
                w = [F(0)]*ncol; w[f] = F(1)
                for i,p in enumerate(pivots):
                    if p is not None and p < ncol: w[p] = -F(str(mat[i,f]))
                ns.append(w)
            # atom images
            v0 = list(vpart)
            for j,c in enumerate(allow):
                if x0[j]:
                    Aj = ACOLS[c]
                    for s in range(NATOM):
                        if Aj[s]: v0[s] += x0[j]*Aj[s]
            NSimg = []
            for w in ns:
                vi = [F(0)]*NATOM
                for j,c in enumerate(allow):
                    if w[j]:
                        Aj = ACOLS[c]
                        for s in range(NATOM):
                            if Aj[s]: vi[s] += w[j]*Aj[s]
                NSimg.append(vi)
            results.append(dict(S=S, v0=v0, NSimg=NSimg, allow=allow, x0=x0))
            mx,dn,nt = height(v0)
            log(f'S={S} CLOSES: nullity={len(ns)} rref-particular max|num|={mx} lcm={dn} terms={nt}')
    json.dump(dict(nsubsets=len(results)),
              open(HERE/'subset_exact_log.json','w'))
    import pickle
    pickle.dump(results, open(HERE/'subset_exact.pkl','wb'))
    log(f'closing subsets: {len(results)}; saved subset_exact.pkl')

if __name__=='__main__':
    main()

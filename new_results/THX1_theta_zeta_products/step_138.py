
import sys, time, pickle, itertools
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, flint, thx1_cert
from math import gcd
import math

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')
def lcm(a, b): return a // gcd(a, b) * b

def zbar(S):
    w = rul1.stuple_to_word(S)
    return {w: Fr((-1)**len(S))}
def zeta_word(n): return {(0,)*(n-1)+(1,): Fr(-1)}
def shpow(w, n):
    out = {w: Fr(1)}
    for _ in range(n-1): out = rul1.sh_mul(out, {w: Fr(1)})
    return out
def li_half(n):
    out = {}
    for c in itertools.product((0, -1), repeat=n-1):
        w = c + (-1,)
        sgn = Fr((-1)**sum(1 for a in c if a == -1))
        out[w] = out.get(w, Fr(0)) + sgn
    return {k: v for k, v in out.items() if v}
LOG2 = {(-1,): Fr(1)}
Z2 = zeta_word(2); Z3 = zeta_word(3); Z5 = zeta_word(5); Z7 = zeta_word(7)
LI4 = li_half(4); LI5 = li_half(5); LI6 = li_half(6); LI7 = li_half(7); LI8 = li_half(8)
X6 = zbar(((5,-1),(1,-1))); X7a = zbar(((5,-1),(1,1),(1,1))); X7b = zbar(((4,-1),(2,1),(1,1)))
X8a = zbar(((6,-1),(1,1),(1,1))); X8b = zbar(((5,-1),(2,1),(1,1)))
X8c = zbar(((5,-1),(1,1),(1,1),(1,1))); X8d = zbar(((4,-1),(2,1),(1,1),(1,1)))
def mono(*factors):
    out = {(): Fr(1)}
    for f in factors: out = rul1.sh_mul(out, f)
    return {k: v for k, v in out.items() if v}
PRIMS = [('log2', 1, LOG2), ('z2', 2, Z2), ('z3', 3, Z3), ('Li4', 4, LI4), ('z5', 5, Z5), ('Li5', 5, LI5),
         ('Li6', 6, LI6), ('X6', 6, X6), ('z7', 7, Z7), ('Li7', 7, LI7), ('X7a', 7, X7a), ('X7b', 7, X7b)]
def monomials_of_weight(W, max_prim_weight):
    prims = [(nm, wt, cb) for nm, wt, cb in PRIMS if wt <= max_prim_weight]
    out = []
    def rec(rem, start, acc):
        if rem == 0:
            out.append(list(acc)); return
        for idx in range(start, len(prims)):
            nm, wt, cb = prims[idx]
            if wt > rem: continue
            acc.append(idx); rec(rem-wt, idx, acc); acc.pop()
    rec(W, 0, [])
    return [("*".join(prims[i][0] for i in inds), mono(*[prims[i][2] for i in inds])) for inds in out]
basis8 = monomials_of_weight(8, 7) + [('Li8', LI8), ('X8a', X8a), ('X8b', X8b), ('X8c', X8c), ('X8d', X8d)]

rels_data = pickle.load(open('rels8_dsdual.pkl','rb'))
rels = [{k: Fr(n, d) for k, (n, d) in r.items()} for fam, tag, r in rels_data]
cols = pickle.load(open('cols8.pkl','rb'))
st = pickle.load(open('struct_w8.pkl','rb'))
B = st['indep_rows']; colpos = st['colpos']; pivcols = st['pivcols']
p0 = 4611686018427387847
r = len(B)
posmap = {c: jj for jj, c in enumerate(pivcols)}
Sq = {}
for i, bi in enumerate(B):
    for k, v in rels[bi].items():
        j = posmap.get(colpos[k])
        if j is not None:
            Sq.setdefault(j, {})[i] = v
L = [1]*r
A = {}
for i in range(r):
    row = Sq.get(i, {})
    for v in row.values():
        L[i] = lcm(L[i], v.denominator)
    A[i] = {j: int(v*L[i]) for j, v in row.items()}
Af = flint.nmod_mat(r, r, p0)
for i, row in A.items():
    for j, v in row.items():
        Af[i, j] = v % p0
Ainv_base = Af.inv()

def I_word_combo(m, j):
    return {(0,)*(m+1) + (-1,)*j: Fr(math.factorial(j)*math.factorial(m))}
def dixon_solve(At, Ainv_t, bb, r, p, nsteps):
    x_hat = [0]*r; pk = 1; bwork = list(bb)
    for step in range(nsteps):
        bv = flint.nmod_mat(r, 1, p)
        for i in range(r):
            bv[i, 0] = bwork[i] % p
        xi = (Ainv_t * bv)
        xil = [int(xi[i, 0]) % p for i in range(r)]
        for i in range(r):
            x_hat[i] += xil[i]*pk
        newb = [0]*r
        for i in range(r):
            s = bwork[i]
            for j, v in At[i].items():
                s -= v*xil[j]
            assert s % p == 0
            newb[i] = s//p
        bwork = newb; pk *= p
    return x_hat, pk

reds8 = pickle.load(open('I_reductions_w89.pkl','rb'))['reds']
for (m, j) in [(0,7), (1,6)]:
    coeffs = reds8[(m,j)]
    V = I_word_combo(m, j)
    for c, (_, cb) in zip(coeffs, basis8):
        if c == 0: continue
        for w, v in cb.items():
            V[w] = V.get(w, Fr(0)) - c*v
    V = {k: v for k, v in V.items() if v}
    b = [Fr(0)]*r
    for k, v in V.items():
        jj = posmap.get(colpos[k])
        if jj is not None:
            b[jj] += v
    Db = 1
    for v in b: Db = lcm(Db, v.denominator)
    At = {i: {jj: v*Db for jj, v in A[i].items()} for i in range(r)}
    bb = [int(b[i]*L[i]*Db) for i in range(r)]
    Ainv_t = Ainv_base * int(pow(Db % p0, -1, p0))
    t=time.time()
    x_hat, pk = dixon_solve(At, Ainv_t, bb, r, p0, 60)
    h = []; fails = 0
    for i in range(r):
        f = thx1_cert.rat_recon(x_hat[i] % pk, pk)
        if f is None: fails += 1
        h.append(f)
    if fails:
        print(f"I_{{{m},{j}}}: recon fails {fails}", flush=True)
        continue
    good, acc = thx1_cert.exact_verify(h, V, rels, B)
    print(f"I_{{{m},{j}}} (w8): {'CERTIFIED' if good else 'FAIL'} ({time.time()-t:.0f}s)", flush=True)
    if good:
        pickle.dump(h, open(f'cert_I_{m}_{j}.pkl','wb'))
        log("I_{%d,%d} reduction CERTIFIED (w8)" % (m, j))

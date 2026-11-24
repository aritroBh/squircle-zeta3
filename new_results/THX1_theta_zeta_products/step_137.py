
import sys, time, pickle, itertools
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, flint, thx1_cert
from math import gcd
import math
import mpmath as mp
mp.mp.dps = 120

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
free_cols = [c for c in range(len(cols)) if c not in posmap]
freepos = {c: jj for jj, c in enumerate(free_cols)}
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
t=time.time()
Ainv = Af.inv()
print("Ainv:", time.time()-t, flush=True)
RBf = flint.nmod_mat(r, 34, p0)
for i, bi in enumerate(B):
    for k, v in rels[bi].items():
        j = freepos.get(colpos[k])
        if j is not None:
            RBf[i, j] = (v.numerator % p0) * pow(v.denominator % p0, -1, p0) % p0
def free_reduce(combo):
    b = flint.nmod_mat(r, 1, p0)
    for k, v in combo.items():
        j = posmap.get(colpos[k])
        if j is not None:
            b[j, 0] = (v.numerator % p0) * pow(v.denominator % p0, -1, p0) % p0
    for i in range(r):
        b[i, 0] = (int(b[i,0]) * (L[i] % p0)) % p0
    x = Ainv * b
    xtRBf = x.transpose() * RBf
    out = [int(xtRBf[0, c]) for c in range(34)]
    for k, v in combo.items():
        j = freepos.get(colpos[k])
        if j is not None:
            out[j] = (out[j] - (v.numerator % p0) * pow(v.denominator % p0, -1, p0) % p0) % p0
    return [(-x) % p0 for x in out]
brows = [free_reduce(combo) for nm, combo in basis8]
Mb = flint.nmod_mat(34, 34, p0)
for i, row in enumerate(brows):
    for c in range(34):
        Mb[c, i] = row[c]
print("basis rank:", Mb.rank(), flush=True)
def I_word_combo(m, j):
    return {(0,)*(m+1) + (-1,)*j: Fr(math.factorial(j)*math.factorial(m))}
def rat_recon(x, P):
    if x == 0: return Fr(0)
    r0, r1 = P, x % P
    t0, t1 = 0, 1
    while r1 * r1 * 2 > P:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        t0, t1 = t1, t0 - q * t1
    if t1 < 0: r1, t1 = -r1, -t1
    if r1 * r1 * 2 > P or t1 * t1 * 2 > P: return None
    return Fr(r1, t1)

l2 = mp.log(2)
z = {n: mp.zeta(n) for n in range(2, 11)}
Li = {n: mp.polylog(n, mp.mpf(1)/2) for n in range(4, 11)}
val = pickle.load(open('basis_vals_89.pkl','rb'))
prim_vals = dict([(nm, v) for nm, wt, v in [('log2',1,l2),('z2',2,z[2]),('z3',3,z[3]),('Li4',4,Li[4]),('z5',5,z[5]),('Li5',5,Li[5]),
    ('Li6',6,Li[6]),('X6',6,val['X6']),('z7',7,z[7]),('Li7',7,Li[7]),('X7a',7,val['X7a']),('X7b',7,val['X7b'])]])
extra_vals = {'Li8': Li[8], 'X8a': val['X8a'], 'X8b': val['X8b'], 'X8c': val['X8c'], 'X8d': val['X8d']}
def name_val(nm):
    if nm in extra_vals: return extra_vals[nm]
    prod = mp.mpf(1)
    for f in nm.split('*'):
        prod *= prim_vals[f]
    return prod
import rul1num
def E1(w, dps=110):
    return rul1num.g_word_over_0m1_num(w, dps=dps)

reds8 = pickle.load(open('I_reductions_w89.pkl','rb'))['reds']
basis_names = [nm for nm, _ in basis8]
for (m, j) in [(0,7), (1,6)]:
    Iv = free_reduce(I_word_combo(m, j))
    bv = flint.nmod_mat(34, 1, p0)
    for c in range(34):
        bv[c, 0] = Iv[c]
    x = Mb.solve(bv)
    coeffs = [rat_recon(int(x[c,0]) % p0, p0) for c in range(34)]
    if any(c is None for c in coeffs):
        print(f"I_{{{m},{j}}}: recon fails", flush=True)
        continue
    Ival = mp.factorial(j)*mp.factorial(m)*E1((0,)*(m+1)+(-1,)*j)
    res = Ival - sum(c*name_val(nm) for c, nm in zip(coeffs, basis_names))
    # compare with PSLQ candidate
    pslq_coeffs = reds8[(m,j)]
    same = all(a == b for a, b in zip(coeffs, pslq_coeffs))
    reds8[(m,j)] = coeffs
    print(f"I_{{{m},{j}}}: |resid| {mp.nstr(abs(res),4)}, matches PSLQ: {same}", flush=True)
pickle.dump({'reds': reds8, 'basis8': basis_names, 'basis9': pickle.load(open('I_reductions_w89.pkl','rb'))['basis9']}, open('I_reductions_w89.pkl','wb'))
log("w8 I_{{0,7}},I_{{1,6}} redone via modular route")

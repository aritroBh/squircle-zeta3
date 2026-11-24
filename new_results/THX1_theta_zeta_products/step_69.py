
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, thx1_cert, thx1_dixon
import sympy as sp
import math, flint
from math import gcd

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

def lcm(a, b): return a // gcd(a, b) * b

def dixon_solve_ainv(A, Ainv, bb, r, p, nsteps):
    x_hat = [0]*r
    pk = 1
    bwork = list(bb)
    for step in range(nsteps):
        bv = flint.nmod_mat(r, 1, p)
        for i in range(r):
            bv[i, 0] = bwork[i] % p
        xi = (Ainv * bv)
        xil = [int(xi[i, 0]) % p for i in range(r)]
        for i in range(r):
            x_hat[i] += xil[i]*pk
        newb = [0]*r
        for i in range(r):
            s = bwork[i]
            for j, v in A[i].items():
                s -= v*xil[j]
            assert s % p == 0, f"residual not divisible at step {step} row {i}"
            newb[i] = s//p
        bwork = newb
        pk *= p
    return x_hat, pk

rels8_data = pickle.load(open('rels8_dsdual.pkl','rb'))
rels = [{k: Fr(n, d) for k, (n, d) in r.items()} for fam, tag, r in rels8_data]
cols8 = pickle.load(open('cols8.pkl','rb'))
st = pickle.load(open('struct_w8.pkl','rb'))
B = st['indep_rows']; colpos = st['colpos']; pivcols = st['pivcols']
r = len(B)
posmap = {c: jj for jj, c in enumerate(pivcols)}
Sq = {}
for i, bi in enumerate(B):
    for k, v in rels[bi].items():
        j = posmap.get(colpos[k])
        if j is not None:
            Sq.setdefault(j, {})[i] = v
# base integer matrix: A_base[i][j] = L_i * Sq[i][j], L_i = lcm of row denominators
L = [1]*r
A = {}
for i in range(r):
    row = Sq.get(i, {})
    for v in row.values():
        L[i] = lcm(L[i], v.denominator)
    A[i] = {j: int(v*L[i]) for j, v in row.items()}
p0 = 4611686018427387847
t=time.time()
Af = flint.nmod_mat(r, r, p0)
for i, row in A.items():
    for j, v in row.items():
        Af[i, j] = v % p0
Ainv_base = Af.inv()
print("Ainv_base:", time.time()-t, flush=True)

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
def V_target(n, m, zb):
    out = Dtilde_combo(n, m)
    coeffs = [Fr(c) for c in Dtarget[f"{n},{m}"]]
    for c, z in zip(coeffs, zb):
        for w, vv in z.items():
            out[w] = out.get(w, Fr(0)) - c*vv
    return {k: v for k, v in out.items() if v}
zb8 = [zprod(8), zprod(3,5), zprod(2,3,3)]

nsteps = 70
for m in range(4):
    V = V_target(8, m, zb8)
    b = [Fr(0)]*r
    for k, v in V.items():
        j = posmap.get(colpos[k])
        if j is not None:
            b[j] += v
    Db = 1
    for v in b: Db = lcm(Db, v.denominator)
    At = {i: {j: v*Db for j, v in A[i].items()} for i in range(r)}
    bb = [int(b[i]*L[i]*Db) for i in range(r)]
    Ainv_t = Ainv_base * int(pow(Db % p0, -1, p0))
    t=time.time()
    x_hat, pk = dixon_solve_ainv(At, Ainv_t, bb, r, p0, nsteps)
    tlift = time.time()-t
    h = []; fails = 0
    for i in range(r):
        f = thx1_cert.rat_recon(x_hat[i] % pk, pk)
        if f is None: fails += 1
        h.append(f)
    print(f"D~_{{8,{m}}}: lift {tlift:.0f}s recon fails {fails}", flush=True)
    if fails: continue
    t=time.time()
    good, acc = thx1_cert.exact_verify(h, V, rels, B)
    print(f"   EXACT VERIFY: {'PASS' if good else 'FAIL'} ({time.time()-t:.0f}s)")
    if good:
        pickle.dump(h, open(f'cert_D8_{m}.pkl','wb'))
        log("D~_{8,%d} PROVED" % m)


import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1
import math, flint
from math import gcd

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')
def lcm(a, b): return a // gcd(a, b) * b

rels9_data = pickle.load(open('rels9_dsdual.pkl','rb'))
rels = [{k: Fr(n, d) for k, (n, d) in r.items()} for fam, tag, r in rels9_data]
cols9 = pickle.load(open('cols9.pkl','rb'))
st = pickle.load(open('struct_w9.pkl','rb'))
B = st['indep_rows']; colpos = st['colpos']; pivcols = st['pivcols']
r = len(B)
posmap = {c: jj for jj, c in enumerate(pivcols)}
Sq = {}
for i, bi in enumerate(B):
    for k, v in rels[bi].items():
        j = posmap.get(colpos[k])
        if j is not None:
            Sq.setdefault(j, {})[i] = v
L = pickle.load(open('w9_L.pkl','rb'))['L']
A = {}
for i in range(r):
    row = Sq.get(i, {})
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
zb9 = [zprod(9), zprod(2,7), zprod(2,2,5), zprod(2,2,2,3), zprod(3,3,3)]

for m in range(4):
    lift = pickle.load(open(f'w9_lift_m{m}.pkl','rb'))
    x_hat, pk, bwork = lift['x_hat'], lift['pk'], lift['bwork']
    # recompute At for this target (need Db)
    V = V_target(9, m, zb9)
    b = [Fr(0)]*r
    for k, v in V.items():
        j = posmap.get(colpos[k])
        if j is not None:
            b[j] += v
    Db = 1
    for v in b: Db = lcm(Db, v.denominator)
    At = {i: {j: v*Db for j, v in A[i].items()} for i in range(r)}
    Ainv_t = Ainv_base * int(pow(Db % p0, -1, p0))
    t=time.time()
    for step in range(70):
        bv = flint.nmod_mat(r, 1, p0)
        for i in range(r):
            bv[i, 0] = bwork[i] % p0
        xi = (Ainv_t * bv)
        xil = [int(xi[i, 0]) % p0 for i in range(r)]
        for i in range(r):
            x_hat[i] += xil[i]*pk
        newb = [0]*r
        for i in range(r):
            s = bwork[i]
            for j, v in At[i].items():
                s -= v*xil[j]
            assert s % p0 == 0
            newb[i] = s//p0
        bwork = newb
        pk *= p0
    pickle.dump({'x_hat': x_hat, 'pk': pk, 'bwork': bwork}, open(f'w9_lift_m{m}.pkl','wb'))
    print(f"m={m}: extended to 160 steps, pk bits {pk.bit_length()}, {time.time()-t:.0f}s", flush=True)
    log("w9 lift m=%d extended to 160 steps" % m)

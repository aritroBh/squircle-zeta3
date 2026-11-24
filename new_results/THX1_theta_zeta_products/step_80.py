
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1
import math, flint
from math import gcd
import numpy as np

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')
def lcm(a, b): return a // gcd(a, b) * b

rels9_data = pickle.load(open('rels9_dsdual.pkl','rb'))
rels = [{k: Fr(n, d) for k, (n, d) in r.items()} for fam, tag, r in rels9_data]
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
arr = np.load('w9_Ainv.npy')
t=time.time()
Ainv_base = flint.nmod_mat(r, r, p0)
for i in range(r):
    for j in range(r):
        Ainv_base[i,j] = int(arr[i,j])
print("Ainv loaded:", time.time()-t, flush=True)

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

def extend_lift(m, target_steps):
    lift = pickle.load(open(f'w9_lift_m{m}.pkl','rb'))
    x_hat, pk, bwork = lift['x_hat'], lift['pk'], lift['bwork']
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
    steps_done = 0
    tmp = pk
    while tmp > 1:
        tmp //= p0; steps_done += 1
    for step in range(steps_done, target_steps):
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
    print(f"m={m}: lift to {target_steps}, pk bits {pk.bit_length()}, {time.time()-t:.0f}s", flush=True)
    log("w9 lift m=%d to step %d" % (m, target_steps))

extend_lift(2, 260)
extend_lift(3, 260)

# sample recon cost
lift = pickle.load(open('w9_lift_m0.pkl','rb'))
x_hat, pk = lift['x_hat'], lift['pk']
def rat_recon_fast(x, P):
    if x == 0: return Fr(0)
    B = P >> 1
    r0, r1 = P, x % P
    t0, t1 = 0, 1
    while True:
        if r1.bit_length()*2 <= B.bit_length()-2 and r1*r1 <= B:
            break
        q, r2 = divmod(r0, r1)
        r0, r1 = r1, r2
        t0, t1 = t1, t0 - q*t1
    if t1 < 0: r1, t1 = -r1, -t1
    if r1*r1*2 > P or t1*t1*2 > P: return None
    return Fr(r1, t1)
t=time.time()
fails = 0
for i in range(0, 8693, 29):
    f = rat_recon_fast(x_hat[i] % pk, pk)
    if f is None: fails += 1
dt = time.time()-t
print(f"sample recon: {fails} fails / 300, {dt:.1f}s -> est {dt*8693/300/60:.1f} min/target")

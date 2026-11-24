
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, thx1_cert
import math

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

rels9_data = pickle.load(open('rels9_dsdual.pkl','rb'))
rels = [{k: Fr(n, d) for k, (n, d) in r.items()} for fam, tag, r in rels9_data]
st = pickle.load(open('struct_w9.pkl','rb'))
B = st['indep_rows']
r = len(B)

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

def rat_recon_fast(x, P):
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

import os
m = 0
lift = pickle.load(open(f'w9_lift_m{m}.pkl','rb'))
x_hat, pk = lift['x_hat'], lift['pk']
t=time.time()
h = []; fails = 0
for i in range(r):
    f = rat_recon_fast(x_hat[i] % pk, pk)
    if f is None: fails += 1
    h.append(f)
    if i % 2000 == 0: print(i, "fails", fails, "t", time.time()-t, flush=True)
print(f"D~_{{9,{m}}}: recon fails {fails}/{r} ({time.time()-t:.0f}s)", flush=True)
if not fails:
    V = V_target(9, m, zb9)
    t=time.time()
    good, acc = thx1_cert.exact_verify(h, V, rels, B)
    print(f"   EXACT VERIFY: {'PASS' if good else 'FAIL'} ({time.time()-t:.0f}s)", flush=True)
    if good:
        pickle.dump(h, open(f'cert_D9_{m}.pkl','wb'))
        log("D~_{9,%d} PROVED" % m)

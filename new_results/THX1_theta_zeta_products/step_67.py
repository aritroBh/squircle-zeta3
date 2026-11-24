
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, thx1_cert, thx1_dixon
import sympy as sp
import math

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

# build w8 relations (exact) and structure
rels8_data = pickle.load(open('rels8_dsdual.pkl','rb'))
rels = [{k: Fr(n, d) for k, (n, d) in r.items()} for fam, tag, r in rels8_data]
cols8 = sorted(rul1.conv_words_of_weight(8))
pickle.dump(cols8, open('cols8.pkl','wb'))
p0 = 4611686018427387847
t=time.time()
st = thx1_cert.compute_structure('w8', rels, cols8, p0, recompute=True)
print(f"w8 structure: rank {st['rank']}, codim {len(cols8)-st['rank']}, {time.time()-t:.0f}s", flush=True)
log("w8 structure: codim %d (DG 34)" % (len(cols8)-st['rank']))

# targets n=8
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

B = st['indep_rows']; colpos = st['colpos']; pivcols = st['pivcols']
r = len(B)
posmap = {c: jj for jj, c in enumerate(pivcols)}
Sq = {}
for i, bi in enumerate(B):
    for k, v in rels[bi].items():
        j = posmap.get(colpos[k])
        if j is not None:
            Sq.setdefault(j, {})[i] = v
t=time.time()
A0, _ = thx1_dixon.prep_system(Sq, [Fr(0)]*r, r)
print("prep A:", time.time()-t)
for m in range(4):
    V = V_target(8, m, zb8)
    b = [Fr(0)]*r
    for k, v in V.items():
        j = posmap.get(colpos[k])
        if j is not None:
            b[j] += v
    A, bb = thx1_dixon.prep_system(Sq, b, r)
    t=time.time()
    x_hat, pk = thx1_dixon.dixon_solve(A, bb, r, p0, target_bits=1100)
    tlift = time.time()-t
    h = []
    fails = 0
    for i in range(r):
        f = thx1_cert.rat_recon(x_hat[i] % pk, pk)
        if f is None: fails += 1
        h.append(f)
    print(f"D~_{{8,{m}}}: lift {tlift:.0f}s pk_bits {pk.bit_length()} recon fails {fails}", flush=True)
    if fails: 
        print("  need more steps")
        continue
    t=time.time()
    good, acc = thx1_cert.exact_verify(h, V, rels, B)
    print(f"   EXACT VERIFY: {'PASS' if good else 'FAIL'} ({time.time()-t:.0f}s)")
    if good:
        pickle.dump(h, open(f'cert_D8_{m}.pkl','wb'))
        log("D~_{8,%d} PROVED" % m)

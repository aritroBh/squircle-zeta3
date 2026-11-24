
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
def V_target(n, m):
    out = Dtilde_combo(n, m)
    coeffs = [Fr(c) for c in Dtarget[f"{n},{m}"]]
    for c, zb in zip(coeffs, [zprod(7), zprod(2,5), zprod(2,2,3)]):
        for w, vv in zb.items():
            out[w] = out.get(w, Fr(0)) - c*vv
    return {k: v for k, v in out.items() if v}

rels_data = pickle.load(open('rels7_dsdual.pkl','rb'))
rels = [{k: Fr(n, d) for k, (n, d) in r.items()} for fam, tag, r in rels_data]
cols7 = sorted(rul1.conv_words_of_weight(7))
st = pickle.load(open('struct_w7.pkl','rb'))
B = st['indep_rows']; colpos = st['colpos']; pivcols = st['pivcols']
r = len(B)
posmap = {c: jj for jj, c in enumerate(pivcols)}
Sq = {}
for i, bi in enumerate(B):
    for k, v in rels[bi].items():
        j = posmap.get(colpos[k])
        if j is not None:
            Sq.setdefault(j, {})[i] = v
V = V_target(7, 0)
b = [Fr(0)]*r
for k, v in V.items():
    j = posmap.get(colpos[k])
    if j is not None:
        b[j] += v

t=time.time()
A, bb = thx1_dixon.prep_system(Sq, b, r)
print("prep:", time.time()-t)
p = 4611686018427387847
t=time.time()
x_hat, pk = thx1_dixon.dixon_solve(A, bb, r, p, target_bits=300)
print("dixon lift:", time.time()-t, "pk bits:", pk.bit_length())

# reconstruct
t=time.time()
h = []
fails = 0
for i in range(r):
    f = thx1_cert.rat_recon(x_hat[i] % pk, pk)
    if f is None: fails += 1
    h.append(f)
print("reconstruct:", time.time()-t, "fails:", fails)

# exact verify
good, acc = thx1_cert.exact_verify(h, V, rels, B)
print("EXACT VERIFY:", "PASS" if good else f"FAIL ({len(acc)})")
if good:
    # compare with saved cert
    cert = pickle.load(open('cert_D7_0.pkl','rb'))
    same = all(h[i] == (sp.Rational(cert[i]) if i in cert else 0) for i in range(r)) if cert else None
    print("matches sympy-solve certificate:", same)
    log("Dixon lifting validated at w7 (D~_{7,0}), exact verify PASS")

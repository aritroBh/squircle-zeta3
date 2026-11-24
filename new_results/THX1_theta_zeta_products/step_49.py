
# fix rref_pivots in thx1_cert.py
src = open('thx1_cert.py').read()
src = src.replace("""def rref_pivots(A):
    R = A.rref()
    piv = []""", """def rref_pivots(A):
    R, rank = A.rref()
    piv = []""")
open('thx1_cert.py','w').write(src)
import importlib, sys
sys.dont_write_bytecode = True
import thx1_cert, rul1, time, pickle
from fractions import Fraction as Fr
import sympy as sp
import math

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

def zeta_word(n):
    return {(0,)*(n-1)+(1,): Fr(-1)}
def zprod(*ns):
    out = {(): Fr(1)}
    for n in ns:
        out = rul1.sh_mul(out, zeta_word(n))
    return {k: v for k, v in out.items() if v}

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
        for w, v in I_word_combo(n-1-j, j).items():
            out[w] = out.get(w, Fr(0)) + c*v
    return {k: v for k, v in out.items() if v}
Dtarget = pickle.load(open('Dtilde_relations.pkl','rb'))
def zbasis_combo(n):
    if n == 2: return [zprod(2)]
    if n == 3: return [zprod(3)]
    if n == 4: return [zprod(4)]
    if n == 5: return [zprod(5), zprod(2,3)]
    if n == 6: return [zprod(6), zprod(3,3)]
    if n == 7: return [zprod(7), zprod(2,5), zprod(2,2,3)]
    if n == 8: return [zprod(8), zprod(3,5), zprod(2,3,3)]
    if n == 9: return [zprod(9), zprod(2,7), zprod(2,2,5), zprod(2,2,2,3), zprod(3,3,3)]
def V_target(n, m):
    out = Dtilde_combo(n, m)
    coeffs = [Fr(c) for c in Dtarget[f"{n},{m}"]]
    for c, zb in zip(coeffs, zbasis_combo(n)):
        for w, v in zb.items():
            out[w] = out.get(w, Fr(0)) - c*v
    return {k: v for k, v in out.items() if v}

primes = [int(sp.prevprime(2**62)), int(sp.prevprime(2**61)), int(sp.prevprime(2**60)), int(sp.prevprime(2**59))]
results = {}
for W in [5, 6, 7]:
    t0 = time.time()
    rels_data = rul1.relations_at_weight(W, families=('DS','DUAL'))
    rels = [r for fam, tag, r in rels_data]
    cols = sorted(rul1.conv_words_of_weight(W))
    st = thx1_cert.compute_structure(f'w{W}', rels, cols, primes[0], recompute=True)
    print(f"w={W}: rels={len(rels)} rank={st['rank']} codim={len(cols)-st['rank']} struct {time.time()-t0:.0f}s", flush=True)
    for m in range(0, W//2):
        V = V_target(W, m)
        t = time.time()
        h, ok = thx1_cert.certify(f'D{W},{m}', V, rels, st, cols, primes[:3])
        print(f"  D~_{{{W},{m}}}: {'EXACT-VERIFIED' if ok else 'FAIL'} ({time.time()-t:.0f}s)", flush=True)
        results[(W,m)] = ok
    log("w=%d done (%.0fs): %s" % (W, time.time()-t0, {k:v for k,v in results.items() if k[0]==W}))
print(results)

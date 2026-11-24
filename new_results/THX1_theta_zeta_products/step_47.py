
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, thx1_cert
import sympy as sp

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

# ---------- building blocks ----------
def zeta_word(n):
    return {(0,)*(n-1)+(1,): Fr(-1)}   # zeta(n) = -G(0^{n-1},1)

def sh_mul(A, B): return rul1.sh_mul(A, B)

def zprod(*ns):
    """zeta(n1)*zeta(n2)*... as word-combo"""
    out = {(): Fr(1)}
    for n in ns:
        out = sh_mul(out, zeta_word(n))
    return {k: v for k, v in out.items() if v}

# ---------- targets: Dtilde combos ----------
e_rm = pickle.load(open('e_rm_ext.pkl','rb'))
def E(r, m): return Fr(e_rm.get(f"{r},{m}", '0'))

def I_word_combo(m, j):
    """I_{m,j} = j! m! G(0^{m+1}(-1)^j)"""
    return {(0,)*(m+1) + (-1,)*j: Fr(m.factorial()*j.factorial() if False else 0)}
import math
def I_word_combo(m, j):
    return {(0,)*(m+1) + (-1,)*j: Fr(math.factorial(j)*math.factorial(m))}

def Dtilde_combo(n, m):
    """2 sum_j e_{m,n-1-j} I_{n-1-j,j}/j! as exact word-combo"""
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
    """Dtilde_{n,m} - (zeta products) as exact word-combo"""
    out = Dtilde_combo(n, m)
    coeffs = [Fr(c) for c in Dtarget[f"{n},{m}"]]
    for c, zb in zip(coeffs, zbasis_combo(n)):
        for w, v in zb.items():
            out[w] = out.get(w, Fr(0)) - c*v
    return {k: v for k, v in out.items() if v}

# ---------- run weights 5,6,7 ----------
primes = [int(sp.prevprime(2**62)), int(sp.prevprime(2**61)), int(sp.prevprime(2**60)), int(sp.prevprime(2**59))]
results = {}
for W in [5, 6, 7]:
    t0 = time.time()
    rels_data = rul1.relations_at_weight(W, families=('DS','DUAL'))
    rels = [r for fam, tag, r in rels_data]
    cols = sorted(rul1.conv_words_of_weight(W))
    st = thx1_cert.compute_structure(f'w{W}', rels, cols, primes[0], recompute=True)
    print(f"w={W}: rels={len(rels)} rank={st['rank']} codim={len(cols)-st['rank']} struct {time.time()-t0:.0f}s", flush=True)
    for n in range(2, W+1):
        for m in range(0, n//2):
            if n != W: continue
            V = V_target(n, m)
            t = time.time()
            h, ok = thx1_cert.certify(f'D{n},{m}', V, rels, st, cols, primes[:3])
            print(f"  D~_{{{n},{m}}}: certificate={'EXACT-VERIFIED' if ok else 'FAIL'} ({time.time()-t:.0f}s)", flush=True)
            results[(n,m)] = ok
    log("w=%d structure+targets done (%.0fs)" % (W, time.time()-t0))
print(results)

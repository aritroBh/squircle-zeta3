
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, thx1_cert
import sympy as sp
import math

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
def zeta_word(n): return {(0,)*(n-1)+(1,): Fr(-1)}
def zprod(*ns):
    out = {(): Fr(1)}
    for n in ns: out = rul1.sh_mul(out, zeta_word(n))
    return {k: v for k, v in out.items() if v}
Dtarget = pickle.load(open('Dtilde_relations.pkl','rb'))
def zbasis_combo(n):
    if n == 7: return [zprod(7), zprod(2,5), zprod(2,2,3)]
def V_target(n, m):
    out = Dtilde_combo(n, m)
    coeffs = [Fr(c) for c in Dtarget[f"{n},{m}"]]
    for c, zb in zip(coeffs, zbasis_combo(n)):
        for w, v in zb.items():
            out[w] = out.get(w, Fr(0)) - c*v
    return {k: v for k, v in out.items() if v}

primes = [int(sp.prevprime(2**62)), int(sp.prevprime(2**61)), int(sp.prevprime(2**60)), int(sp.prevprime(2**59)), int(sp.prevprime(2**58))]
W = 7
rels_data = rul1.relations_at_weight(W, families=('DS','DUAL'))
rels = [r for fam, tag, r in rels_data]
cols = sorted(rul1.conv_words_of_weight(W))
st = pickle.load(open('struct_w7.pkl','rb'))
V = V_target(7, 0)
print("target size:", len(V))

# solve at 5 primes, check reconstruction stats
hmods = []
for p in primes:
    hs = thx1_cert.solve_h_mod_batch([V], rels, st, cols, p)
    assert thx1_cert.check_h_mod(hs[0], V, rels, st, cols, p)
    hmods.append(hs[0])
print("mod-p solves + checks OK at 5 primes")
r = len(st['indep_rows'])
# how many entries reconstruct with k primes?
for k in range(2, 6):
    ok = 0; fail = 0
    for i in range(r):
        x, P = thx1_cert.crt([(hmods[j][i], primes[j]) for j in range(k)])
        f = thx1_cert.rat_recon(x, P)
        if f is None: fail += 1
        else: ok += 1
    print(f"k={k} primes: reconstructed {ok}, failed {fail}")

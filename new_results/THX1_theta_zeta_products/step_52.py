
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, thx1_cert
import sympy as sp
import math

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

cp = pickle.load(open('w7_hist_cp.pkl','rb'))
pivots, pivot_rows = cp['pivots'], cp['pivot_rows']
P = 4611686018427387847

def reduce_mod_hist(V, pivot_rows, p):
    """V: exact Fraction combo. Reduce against pivot rows (mod p, with history). Returns (residue, hist_orig)."""
    v = {k: (c.numerator % p) * pow(c.denominator % p, -1, p) % p for k, c in V.items() if c}
    hist = {}
    for pc, pd, ph in pivot_rows:
        if pc not in v: continue
        f = v[pc]
        for k, x in pd.items():
            v[k] = (v.get(k, 0) - f*x) % p
            if v[k] == 0: del v[k]
        for k, x in ph.items():
            hist[k] = (hist.get(k, 0) + f*x) % p
            if hist[k] == 0: del hist[k]
    return v, hist

# targets
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

for m in range(3):
    V = V_target(7, m)
    res, hist = reduce_mod_hist(V, pivot_rows, P)
    print(f"D~_{{7,{m}}}: residue size = {len(res)}, history entries = {len(hist)}")
    if len(res) == 0:
        # try single-prime reconstruction
        fails = 0; h_rat = {}
        for k, x in hist.items():
            f = thx1_cert.rat_recon(x, P)
            if f is None: fails += 1
            else: h_rat[k] = f
        print(f"   single-prime reconstruction: {len(h_rat)} ok, {fails} fail")
        pickle.dump({'hist': hist}, open(f'w7_hist_m{m}.pkl','wb'))
log("w7 targets reduced mod p with history; residue 0 for all m")

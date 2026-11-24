
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, thx1_cert
import math

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

primes = [4611686018427387847, 2305843009213693951, 1152921504606846883,
          576460752303423433, 288230376151711717, 144115188075855859]
cps = [pickle.load(open('w7_hist_cp.pkl','rb'))] + [pickle.load(open(f'w7_hist_p{i}.pkl','rb')) for i in range(2,7)]

def reduce_mod_hist(V, pivot_rows, p):
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

for m in range(3):
    V = V_target(7, m)
    hists = []
    for i, cp in enumerate(cps):
        res, hist = reduce_mod_hist(V, cp['pivot_rows'], primes[i])
        assert len(res) == 0
        hists.append(hist)
    # per-entry: reconstruct using k primes, k = 2..6; take first that verifies mod next prime
    h_rat = {}
    worst = 0
    for key in hists[0]:
        got = None
        for k in range(2, 7):
            x, P = thx1_cert.crt([(hists[j][key], primes[j]) for j in range(k)])
            f = thx1_cert.rat_recon(x, P)
            if f is not None:
                got = f; break
        if got is None:
            h_rat = None; break
        h_rat[key] = got
        worst = max(worst, max(abs(got.numerator), got.denominator).bit_length())
    if h_rat is None:
        print(f"D~_{{7,{m}}}: STILL FAILING entries after 6 primes")
        continue
    t=time.time()
    good, acc = thx1_cert.exact_verify([h_rat.get(i, Fr(0)) for i in range(len(rels))], V, rels, list(range(len(rels))))
    print(f"D~_{{7,{m}}}: reconstructed {len(h_rat)} entries (max height {worst} bits), EXACT VERIFY: {'PASS' if good else 'FAIL'} ({time.time()-t:.0f}s)")
    if good:
        pickle.dump(h_rat, open(f'cert_D7_{m}.pkl','wb'))
        log("D~_{7,%d} PROVED with exact certificate (%d instances, max height %d bits)" % (m, len(h_rat), worst))

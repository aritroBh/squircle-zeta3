
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, thx1_cert, itertools, math
import sympy as sp

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

def zbar(S):
    w = rul1.stuple_to_word(S)
    return {w: Fr((-1)**len(S))}
def zeta_word(n): return {(0,)*(n-1)+(1,): Fr(-1)}
def shpow(w, n):
    out = {w: Fr(1)}
    for _ in range(n-1): out = rul1.sh_mul(out, {w: Fr(1)})
    return out
def li_half(n):
    out = {}
    for c in itertools.product((0, -1), repeat=n-1):
        w = c + (-1,)
        sgn = Fr((-1)**sum(1 for a in c if a == -1))
        out[w] = out.get(w, Fr(0)) + sgn
    return {k: v for k, v in out.items() if v}
LOG2 = {(-1,): Fr(1)}
Z2 = zeta_word(2); Z3 = zeta_word(3); Z5 = zeta_word(5); Z7 = zeta_word(7); Z9 = zeta_word(9)
LI4 = li_half(4); LI5 = li_half(5); LI6 = li_half(6); LI7 = li_half(7)
def mono(*factors):
    out = {(): Fr(1)}
    for f in factors: out = rul1.sh_mul(out, f)
    return {k: v for k, v in out.items() if v}
X6 = zbar(((5,-1),(1,-1)))
X7a = zbar(((5,-1),(1,1),(1,1)))
X7b = zbar(((4,-1),(2,1),(1,1)))

basis_combos = {
 6: [('log2^6', shpow((-1,),6)), ('z2*log2^4', mono(Z2, shpow((-1,),4))), ('z2^2*log2^2', mono(Z2,Z2,shpow((-1,),2))),
     ('z2^3', mono(Z2,Z2,Z2)), ('z3*log2^3', mono(Z3, shpow((-1,),3))), ('z3*z2*log2', mono(Z3,Z2,LOG2)),
     ('z3^2', mono(Z3,Z3)), ('z5*log2', mono(Z5,LOG2)), ('Li4*log2^2', mono(LI4, shpow((-1,),2))),
     ('Li4*z2', mono(LI4,Z2)), ('Li5*log2', mono(LI5,LOG2)), ('Li6', LI6), ('X6', X6)],
 7: [('log2^7', shpow((-1,),7)), ('z2*log2^5', mono(Z2, shpow((-1,),5))), ('z2^2*log2^3', mono(Z2,Z2,shpow((-1,),3))),
     ('z2^3*log2', mono(Z2,Z2,Z2,LOG2)), ('z3*log2^4', mono(Z3, shpow((-1,),4))), ('z3*z2*log2^2', mono(Z3,Z2,shpow((-1,),2))),
     ('z3*z2^2', mono(Z3,Z2,Z2)), ('z3^2*log2', mono(Z3,Z3,LOG2)), ('z5*log2^2', mono(Z5, shpow((-1,),2))),
     ('z5*z2', mono(Z5,Z2)), ('Li4*log2^3', mono(LI4, shpow((-1,),3))), ('Li4*z2*log2', mono(LI4,Z2,LOG2)),
     ('Li4*z3', mono(LI4,Z3)), ('Li5*log2^2', mono(LI5, shpow((-1,),2))), ('Li5*z2', mono(LI5,Z2)),
     ('Li6*log2', mono(LI6,LOG2)), ('X6*log2', mono(X6,LOG2)), ('z7', Z7), ('Li7', LI7), ('X7a', X7a), ('X7b', X7b)],
}

reds = pickle.load(open('I_reductions_w2_7.pkl','rb'))
primes = [int(sp.prevprime(2**62)), int(sp.prevprime(2**61)), int(sp.prevprime(2**60))]

def I_word_combo(m, j):
    return {(0,)*(m+1) + (-1,)*j: Fr(math.factorial(j)*math.factorial(m))}

certs_I = {}
for W in [6, 7]:
    rels_data = pickle.load(open(f'rels{W}_dsdual.pkl','rb'))
    rels = [{k: Fr(n, d) for k, (n, d) in r.items()} for fam, tag, r in rels_data]
    cols = pickle.load(open(f'cols{W}.pkl','rb'))
    st = pickle.load(open(f'struct_w{W}.pkl','rb'))
    bcombos = basis_combos[W]
    for m in range(W-1):
        j = W-1-m
        coeffs = reds[(m,j)]
        V = I_word_combo(m, j)
        for c, (_, cb) in zip(coeffs, bcombos):
            if c == 0: continue
            for w, v in cb.items():
                V[w] = V.get(w, Fr(0)) - c*v
        V = {k: v for k, v in V.items() if v}
        t=time.time()
        h, ok = thx1_cert.certify(f'I{m},{j}', V, rels, st, cols, primes)
        print(f"I_{{{m},{j}}} (w{W}): {'CERTIFIED' if ok else 'FAIL'} ({time.time()-t:.0f}s)", flush=True)
        if ok:
            certs_I[(m,j)] = h
            pickle.dump(h, open(f'cert_I_{m}_{j}.pkl','wb'))
            log("I_{%d,%d} reduction CERTIFIED" % (m, j))

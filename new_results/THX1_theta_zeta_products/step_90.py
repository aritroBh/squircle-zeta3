
import sys, time, pickle, os
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, flint, itertools

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

# build & save rels for w=5,6 if missing
for W in [5, 6]:
    fn = f'rels{W}_dsdual.pkl'
    if not os.path.exists(fn):
        rels_data = rul1.relations_at_weight(W, families=('DS','DUAL'))
        data = [(fam, tag, {k: (v.numerator, v.denominator) for k, v in r.items()}) for fam, tag, r in rels_data]
        pickle.dump(data, open(fn,'wb'))
        cols = sorted(rul1.conv_words_of_weight(W))
        pickle.dump(cols, open(f'cols{W}.pkl','wb'))
        print(f"saved {fn} ({len(data)} rels)")

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
LI4 = li_half(4); LI5 = li_half(5); LI6 = li_half(6); LI7 = li_half(7); LI8 = li_half(8); LI9 = li_half(9)
def mono(*factors):
    out = {(): Fr(1)}
    for f in factors: out = rul1.sh_mul(out, f)
    return {k: v for k, v in out.items() if v}

def check_basis(W, named):
    rels_data = pickle.load(open(f'rels{W}_dsdual.pkl','rb'))
    rels = [{k: Fr(n, d) for k, (n, d) in r.items()} for fam, tag, r in rels_data]
    cols = pickle.load(open(f'cols{W}.pkl','rb'))
    p0 = 4611686018427387847
    colpos = {c:i for i,c in enumerate(cols)}
    A = flint.nmod_mat(len(rels)+len(named), len(cols), p0)
    for i, r in enumerate(rels):
        for k, v in r.items():
            A[i, colpos[k]] = (v.numerator % p0) * pow(v.denominator % p0, -1, p0) % p0
    for j, (nm, combo) in enumerate(named):
        for k, v in combo.items():
            A[len(rels)+j, colpos[k]] = (v.numerator % p0) * pow(v.denominator % p0, -1, p0) % p0
    return A.rank()

classical6 = [('log2^6', shpow((-1,),6)), ('z2 log2^4', mono(Z2, shpow((-1,),4))), ('z2^2 log2^2', mono(Z2,Z2,shpow((-1,),2))),
              ('z2^3', mono(Z2,Z2,Z2)), ('z3 log2^3', mono(Z3, shpow((-1,),3))), ('z3 z2 log2', mono(Z3,Z2,LOG2)),
              ('z3^2', mono(Z3,Z3)), ('z5 log2', mono(Z5,LOG2)), ('Li4 log2^2', mono(LI4, shpow((-1,),2))),
              ('Li4 z2', mono(LI4,Z2)), ('Li5 log2', mono(LI5,LOG2)), ('Li6', LI6)]
for X6name, X6 in [('zbar(5,1)', zbar(((5,-1),(1,1)))), ('zbar(5,-1)', zbar(((5,-1),(1,-1))))]:
    rk = check_basis(6, classical6 + [(X6name, X6)])
    print(f"w6 with {X6name}: rank={rk} (need 324) {'VALID' if rk==324 else 'invalid'}", flush=True)
log("w6 basis candidates checked")

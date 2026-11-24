
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, flint, itertools

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

classical7 = [('log2^7', shpow((-1,),7)), ('z2 log2^5', mono(Z2, shpow((-1,),5))), ('z2^2 log2^3', mono(Z2,Z2,shpow((-1,),3))),
              ('z2^3 log2', mono(Z2,Z2,Z2,LOG2)), ('z3 log2^4', mono(Z3, shpow((-1,),4))), ('z3 z2 log2^2', mono(Z3,Z2,shpow((-1,),2))),
              ('z3 z2^2', mono(Z3,Z2,Z2)), ('z3^2 log2', mono(Z3,Z3,LOG2)), ('z5 log2^2', mono(Z5, shpow((-1,),2))),
              ('z5 z2', mono(Z5,Z2)), ('Li4 log2^3', mono(LI4, shpow((-1,),3))), ('Li4 z2 log2', mono(LI4,Z2,LOG2)),
              ('Li4 z3', mono(LI4,Z3)), ('Li5 log2^2', mono(LI5, shpow((-1,),2))), ('Li5 z2', mono(LI5,Z2)),
              ('Li6 log2', mono(LI6,LOG2)), ('X6 log2', mono(X6,LOG2)), ('z7', Z7), ('Li7', LI7)]
cands = [('zbar(6,-1)', zbar(((6,-1),(1,-1)))), ('zbar(5,-2)', zbar(((5,-1),(2,-1)))), ('zbar(6,1)', zbar(((6,-1),(1,1)))),
         ('zbar(5,2)', zbar(((5,-1),(2,1)))), ('zbar(5,-1,-1)', zbar(((5,-1),(1,-1),(1,-1)))), ('zbar(5,1,1)', zbar(((5,-1),(1,1),(1,1)))),
         ('zbar(4,-3)', zbar(((4,-1),(3,-1)))), ('zbar(3,-4)', zbar(((3,-1),(4,-1))))]
import itertools as it
found = None
for (n1,c1),(n2,c2) in it.combinations(cands, 2):
    rk = check_basis(7, classical7 + [(n1,c1),(n2,c2)])
    if rk == 972:
        print(f"w7 VALID with {n1}, {n2}")
        found = (n1,n2)
        break
if not found:
    print("no pair worked; need deeper candidates")
log("w7 basis: %s" % (found,))

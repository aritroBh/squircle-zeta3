
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
LI4 = li_half(4); LI5 = li_half(5); LI6 = li_half(6); LI7 = li_half(7); LI8 = li_half(8); LI9 = li_half(9)
X6 = zbar(((5,-1),(1,-1))); X7a = zbar(((5,-1),(1,1),(1,1))); X7b = zbar(((4,-1),(2,1),(1,1)))
X8a = zbar(((6,-1),(1,1),(1,1))); X8b = zbar(((5,-1),(2,1),(1,1)))
X8c = zbar(((5,-1),(1,1),(1,1),(1,1))); X8d = zbar(((4,-1),(2,1),(1,1),(1,1)))
def mono(*factors):
    out = {(): Fr(1)}
    for f in factors: out = rul1.sh_mul(out, f)
    return {k: v for k, v in out.items() if v}
PRIMS = [('log2', 1, LOG2), ('z2', 2, Z2), ('z3', 3, Z3), ('Li4', 4, LI4), ('z5', 5, Z5), ('Li5', 5, LI5),
         ('Li6', 6, LI6), ('X6', 6, X6), ('z7', 7, Z7), ('Li7', 7, LI7), ('X7a', 7, X7a), ('X7b', 7, X7b),
         ('Li8', 8, LI8), ('X8a', 8, X8a), ('X8b', 8, X8b), ('X8c', 8, X8c), ('X8d', 8, X8d)]
def monomials_of_weight(W, max_prim_weight):
    prims = [(nm, wt, cb) for nm, wt, cb in PRIMS if wt <= max_prim_weight]
    out = []
    def rec(rem, start, acc):
        if rem == 0:
            out.append(list(acc)); return
        for idx in range(start, len(prims)):
            nm, wt, cb = prims[idx]
            if wt > rem: continue
            acc.append(idx); rec(rem-wt, idx, acc); acc.pop()
    rec(W, 0, [])
    return [("*".join(prims[i][0] for i in inds), mono(*[prims[i][2] for i in inds])) for inds in out]

classical9 = monomials_of_weight(9, 8) + [('z9', Z9), ('Li9', LI9)]
# candidate final extras: 5 depth3/4 + X9f
cand_extras = [('zbar(7,1,1)', zbar(((7,-1),(1,1),(1,1)))),
               ('zbar(6,2,1)', zbar(((6,-1),(2,1),(1,1)))),
               ('zbar(5,3,1)', zbar(((5,-1),(3,1),(1,1)))),
               ('zbar(5,2,2)', zbar(((5,-1),(2,1),(2,1)))),
               ('zbar(4,4,1)', zbar(((4,-1),(4,1),(1,1)))),
               ('X9f', {(0,1,1,1,1,1,1,0,-1): Fr(1)})]
rels_data = pickle.load(open('rels9_dsdual.pkl','rb'))
rels = [{k: Fr(n, d) for k, (n, d) in r.items()} for fam, tag, r in rels_data]
cols9 = pickle.load(open('cols9.pkl','rb'))
p0 = 4611686018427387847
colpos = {c:i for i,c in enumerate(cols9)}
named = classical9 + cand_extras
t=time.time()
A = flint.nmod_mat(len(rels)+len(named), len(cols9), p0)
for i, r in enumerate(rels):
    for k, v in r.items():
        A[i, colpos[k]] = (v.numerator % p0) * pow(v.denominator % p0, -1, p0) % p0
for j, (nm, combo) in enumerate(named):
    for k, v in combo.items():
        A[len(rels)+j, colpos[k]] = (v.numerator % p0) * pow(v.denominator % p0, -1, p0) % p0
rk = A.rank()
print("final w9 basis rank:", rk, "(need 8748)", "time", time.time()-t, flush=True)
log("w9 final basis rank %d" % rk)

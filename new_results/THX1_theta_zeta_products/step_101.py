
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, thx1_cert, thx1_dixon, itertools, math, flint
from math import gcd

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')
def lcm(a, b): return a // gcd(a, b) * b

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
Z2 = zeta_word(2); Z3 = zeta_word(3); Z5 = zeta_word(5); Z7 = zeta_word(7)
LI4 = li_half(4); LI5 = li_half(5); LI6 = li_half(6); LI7 = li_half(7)
def mono(*factors):
    out = {(): Fr(1)}
    for f in factors: out = rul1.sh_mul(out, f)
    return {k: v for k, v in out.items() if v}
X6 = zbar(((5,-1),(1,-1)))
X7a = zbar(((5,-1),(1,1),(1,1)))
X7b = zbar(((4,-1),(2,1),(1,1)))
bcombos7 = [('log2^7', shpow((-1,),7)), ('z2*log2^5', mono(Z2, shpow((-1,),5))), ('z2^2*log2^3', mono(Z2,Z2,shpow((-1,),3))),
     ('z2^3*log2', mono(Z2,Z2,Z2,LOG2)), ('z3*log2^4', mono(Z3, shpow((-1,),4))), ('z3*z2*log2^2', mono(Z3,Z2,shpow((-1,),2))),
     ('z3*z2^2', mono(Z3,Z2,Z2)), ('z3^2*log2', mono(Z3,Z3,LOG2)), ('z5*log2^2', mono(Z5, shpow((-1,),2))),
     ('z5*z2', mono(Z5,Z2)), ('Li4*log2^3', mono(LI4, shpow((-1,),3))), ('Li4*z2*log2', mono(LI4,Z2,LOG2)),
     ('Li4*z3', mono(LI4,Z3)), ('Li5*log2^2', mono(LI5, shpow((-1,),2))), ('Li5*z2', mono(LI5,Z2)),
     ('Li6*log2', mono(LI6,LOG2)), ('X6*log2', mono(X6,LOG2)), ('z7', Z7), ('Li7', LI7), ('X7a', X7a), ('X7b', X7b)]

rels_data = pickle.load(open('rels7_dsdual.pkl','rb'))
rels = [{k: Fr(n, d) for k, (n, d) in r.items()} for fam, tag, r in rels_data]
cols = pickle.load(open('cols7.pkl','rb'))
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
L = [1]*r
A = {}
for i in range(r):
    row = Sq.get(i, {})
    for v in row.values():
        L[i] = lcm(L[i], v.denominator)
    A[i] = {j: int(v*L[i]) for j, v in row.items()}
p0 = 4611686018427387847
Af = flint.nmod_mat(r, r, p0)
for i, row in A.items():
    for j, v in row.items():
        Af[i, j] = v % p0
Ainv_base = Af.inv()
print("Ainv ready", flush=True)

def I_word_combo(m, j):
    return {(0,)*(m+1) + (-1,)*j: Fr(math.factorial(j)*math.factorial(m))}
reds = pickle.load(open('I_reductions_w2_7.pkl','rb'))

def dixon_solve_ainv(At, Ainv_t, bb, r, p, nsteps):
    x_hat = [0]*r; pk = 1; bwork = list(bb)
    for step in range(nsteps):
        bv = flint.nmod_mat(r, 1, p)
        for i in range(r):
            bv[i, 0] = bwork[i] % p
        xi = (Ainv_t * bv)
        xil = [int(xi[i, 0]) % p for i in range(r)]
        for i in range(r):
            x_hat[i] += xil[i]*pk
        newb = [0]*r
        for i in range(r):
            s = bwork[i]
            for j, v in At[i].items():
                s -= v*xil[j]
            assert s % p == 0
            newb[i] = s//p
        bwork = newb; pk *= p
    return x_hat, pk

for m in range(6):
    j = 6-m
    coeffs = reds[(m,j)]
    V = I_word_combo(m, j)
    for c, (_, cb) in zip(coeffs, bcombos7):
        if c == 0: continue
        for w, v in cb.items():
            V[w] = V.get(w, Fr(0)) - c*v
    V = {k: v for k, v in V.items() if v}
    b = [Fr(0)]*r
    for k, v in V.items():
        jj = posmap.get(colpos[k])
        if jj is not None:
            b[jj] += v
    Db = 1
    for v in b: Db = lcm(Db, v.denominator)
    At = {i: {jj: v*Db for jj, v in A[i].items()} for i in range(r)}
    bb = [int(b[i]*L[i]*Db) for i in range(r)]
    Ainv_t = Ainv_base * int(pow(Db % p0, -1, p0))
    x_hat, pk = dixon_solve_ainv(At, Ainv_t, bb, r, p0, 12)
    h = []; fails = 0
    for i in range(r):
        f = thx1_cert.rat_recon(x_hat[i] % pk, pk)
        if f is None: fails += 1
        h.append(f)
    if fails:
        print(f"I_{{{m},{j}}}: recon fails {fails}")
        continue
    good, acc = thx1_cert.exact_verify(h, V, rels, B)
    print(f"I_{{{m},{j}}} (w7): {'CERTIFIED' if good else 'FAIL'}", flush=True)
    if good:
        pickle.dump(h, open(f'cert_I_{m}_{j}.pkl','wb'))
        log("I_{%d,%d} reduction CERTIFIED (w7)" % (m, j))

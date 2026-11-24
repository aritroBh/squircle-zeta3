
import sys, time, pickle
sys.dont_write_bytecode = True
import rul1, flint, sympy as sp

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

# primitive counts from Hilbert series d = 1,1,2,3,5,8,13,21,34,55
d = [1,1,2,3,5,8,13,21,34,55]
import sympy
# log H = sum a_k t^k ; k a_k = sum_{w|k} w p_w ; mobius invert
mu = {n: int(sympy.mobius(n)) for n in range(1, 12)}
a = [0]*11
H = d + [0]*(11-len(d))
# compute log H coefficients
logc = [0]*11
for k in range(1, 11):
    s = sum(-logc[j]*H[k-j]* (k-j)/k for j in range(1, k))
    logc[k] = H[k] - s
# simpler: a_k from derivative: (log H)' = H'/H
# use sympy
t = sp.symbols('t')
Hs = sum(d[i]*t**i for i in range(len(d)))
logH = sp.series(sp.log(Hs), t, 0, 11).removeO().expand()
p = []
for w in range(1, 10):
    ak = sp.expand(logH).coeff(t, w)
    # w p_w = sum_{dv|w} mu(d) (w/d) a_{w/d}... use mobius on k*a_k
    ka = {k: sp.expand(logH).coeff(t, k)*k for k in range(1, 10)}
    wpw = sum(mu[dv] * ka[w//dv] for dv in range(1, w+1) if w % dv == 0)
    p.append(int(sp.simplify(wpw/w)))
print("primitive counts p_w (w=1..9):", p)
log("level-2 primitive counts w=1..9: %s" % p)

# build w9 relations and matrix, flint rank
t0=time.time()
rels9 = rul1.relations_at_weight(9, families=('DS','DUAL'))
print("build w9 rels: %.0fs, %d rels"%(time.time()-t0, len(rels9)), flush=True)
data = [(fam, tag, {k: (v.numerator, v.denominator) for k, v in r.items()}) for fam, tag, r in rels9]
pickle.dump(data, open('rels9_dsdual.pkl','wb'))
cols9 = sorted(rul1.conv_words_of_weight(9))
colpos9 = {c:i for i,c in enumerate(cols9)}
pickle.dump(cols9, open('cols9.pkl','wb'))
P = int(sp.prevprime(2**62))
t = time.time()
A = flint.nmod_mat(len(data), len(cols9), P)
for i, (fam, tag, r) in enumerate(data):
    for k, (n, dd) in r.items():
        A[i, colpos9[k]] = (n * pow(dd, -1, P)) % P
print("build matrix %.0fs"%(time.time()-t), flush=True)
t = time.time()
r = A.rank()
print("flint rank w9: %d in %.0fs (codim %d, DG 55)" % (r, time.time()-t, len(cols9)-r), flush=True)
log("w9 DS+DUAL codim = %d (DG 55) via flint rank, %.0fs" % (len(cols9)-r, time.time()-t))

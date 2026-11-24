
import pickle, sympy as sp, time
import mpmath as mp
t0 = time.time()
D1 = '../DOD1_dodecagon/'
r3 = sp.sqrt(3)
sad = pickle.load(open(D1+'S_atom_decomp.pkl','rb'))
coefs = sad['coefs']
def to_pq(expr):
    e = sp.radsimp(sp.expand(expr))
    e = sp.expand(e)
    q = e.coeff(r3)
    p = sp.expand(e - q*r3)
    return sp.Rational(p), sp.Rational(q)
out = {}
mp.mp.dps = 120
for key in coefs:
    if key[0] not in ('Lam','T3a'): continue
    p, q = to_pq(coefs[key])
    out[key] = (p, q)
    # high-precision numeric verification of the reconstruction
    v_exact = mp.mpf(sp.N(p,120)) + mp.mpf(sp.N(q,120))*mp.sqrt(3)
    v_orig = mp.mpf(sp.N(coefs[key],120))
    d = abs(v_exact - v_orig)
    print(key, "= (%s)+(%s)sqrt3  [check %.1e]" % (p, q, float(d)))
pickle.dump(out, open('atom_coefs_pq.pkl','wb'))
print("saved atom_coefs_pq.pkl (%d entries), elapsed %ds" % (len(out), time.time()-t0))


import pickle, sympy as sp, time
t0 = time.time()
D1 = '../DOD1_dodecagon/'
r3 = sp.sqrt(3)
sad = pickle.load(open(D1+'S_atom_decomp.pkl','rb'))
coefs = sad['coefs']

def to_pq(expr):
    # expr in Q(sqrt3); return (p,q) rationals with expr = p + q*sqrt(3)
    e = sp.radsimp(sp.expand(expr))
    e = sp.expand(e)
    q = sp.expand(e).coeff(r3)
    p = sp.expand(e - q*r3)
    assert q is not None
    p = sp.Rational(p); q = sp.Rational(q)
    return p, q

out = {}
worst = 0
for key in coefs:
    p, q = to_pq(coefs[key])
    # numeric check
    import mpmath as mp
    mp.mp.dps = 40
    v_exact = float(p) + float(q)*mp.sqrt(3)
    v_orig = sp.N(coefs[key], 40)
    d = abs(v_exact - v_orig)
    worst = max(worst, d)
    out[key] = (p, q)
    print(key, "= (%s) + (%s)*sqrt3   [check %.2e]" % (p, q, d))
print("worst reconstruction error:", worst)
pickle.dump(out, open('atom_coefs_pq.pkl','wb'))
print("saved atom_coefs_pq.pkl, elapsed", time.time()-t0)

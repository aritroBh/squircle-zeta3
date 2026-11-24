
import pickle
import sympy as sp
t = sp.Symbol('t'); s5 = sp.sqrt(5)
mons = pickle.load(open("E5_dec_poly.pkl","rb"))

def qq5_extract(cf):
    cf = sp.radsimp(sp.expand(cf))
    cf = sp.expand(cf)
    A = cf.coeff(s5, 0); B = cf.coeff(s5, 1)
    assert sp.expand(cf - (A + B*s5)) == 0, ("extract fail", cf)
    return sp.Rational(A), sp.Rational(B)

parsed = []
for mon, c in mons:
    cc = sp.sympify(c, locals={'t': t})
    p = sp.Poly(sp.expand(cc), t)
    out = {}
    for (n,), cf in p.terms():
        out[n] = qq5_extract(cf)
    parsed.append((mon, out))
print("parsed %d monomials OK" % len(parsed))
pickle.dump(parsed, open("parsed_mons_dec.pkl","wb"))
# sanity: all coefficients rational (p,q)
maxdeg = max(n for mon,cu in parsed for n in cu)
print("max t-degree:", maxdeg)
need5 = set()
for mon, cu in parsed:
    if sum(mon)==5:
        a,b,c,d,e,f = mon
        for n in cu: need5.add((n,a,b,c,d))
print("weight-5 integrals needed:", len(need5))
needall = set()
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    for n in cu: needall.add((n,a,b,c,d))
print("all integrals needed:", len(needall))
print(parsed[0])

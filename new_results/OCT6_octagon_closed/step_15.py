import pickle, sympy as sp, mpmath as mp
from o6util import load
mons = load('E3_monomials.pkl')['mons']
t, r2 = sp.symbols('t'), sp.sqrt(2)
v = sp.sqrt(10-7*r2); u = sp.sqrt(2+r2)

def to_uq(coeff):
    # coeff = g/v with g in Q(r2); coeff = u*g*(2+r2)/2
    g = sp.radsimp(sp.simplify(coeff*v))
    q = sp.expand(g*(2+r2)/2)
    return sp.simplify(q)

simple = []
for key, coeff in mons:
    tp = sp.Poly(sp.expand(coeff), t)
    terms = [(n[0], to_uq(c)) for n, c in tp.terms()]
    simple.append((key, terms))

# verify: coeff == u*sum(q t^n) numerically for all 83 at t=0.31
mp.mp.dps=40
tt = mp.mpf('0.31'); uu = mp.sqrt(2+mp.sqrt(2))
maxd = mp.mpf(0)
for (key, coeff), (key2, terms) in zip(mons, simple):
    c1 = mp.mpf(sp.N(coeff.subs(t, tt), 50))
    c2 = uu*sum(mp.mpf(sp.N(q.subs(t,tt),50)) if sp.N(q,30)!=0 else 0 for n,q in terms)*tt**0  # q already scalar per term
    c2 = uu*sum(mp.mpf(sp.N(q,50))*tt**n for n,q in terms)
    maxd = max(maxd, abs(c1-c2))
print("max |coeff - u*q(t)| over 83 mons:", mp.nstr(maxd,3))

# print the 83 monomials: key=(a,b,c1,c2,e1,e2) exponents of (log t, log(1-t), W1, W2, c1, c2)
print("\n83 monomials  E3 = u * Σ q(t)·X^a Y^b W1^c1 W2^c2 ·c1^e1 c2^e2,  u=sqrt(2+√2), c1=-l2/2-lu, c2=-lu")
for key, terms in simple:
    qs = ' + '.join(f"({sp.nsimplify(q)})·t^{n}" if n else f"({sp.nsimplify(q)})" for n,q in sorted(terms, reverse=True))
    print(key, ':', qs)
pickle.dump(simple, open('mons_simple.pkl','wb'))

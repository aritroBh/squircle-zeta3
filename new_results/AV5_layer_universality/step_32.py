import mpmath, pickle, itertools, json
from mpmath import mp, mpf, sqrt, zeta, log, polylog, nstr
mp.dps = 90
gen = [(1,'L',log(2)),(2,'Z2',zeta(2)),(3,'Z3',zeta(3)),(4,'Li4',polylog(4,mpf('0.5'))),
       (5,'Z5',zeta(5)),(5,'Li5',polylog(5,mpf('0.5'))),(6,'Li6',polylog(6,mpf('0.5'))),
       (7,'Z7',zeta(7)),(7,'Li7',polylog(7,mpf('0.5')))]
def monomials(w):
    """all monomials of exact weight w"""
    out = []
    def rec(rem, start, factors, val):
        if rem==0:
            out.append((factors, val)); return
        for i in range(start, len(gen)):
            gw, gn, gv = gen[i]
            if gw <= rem:
                rec(rem-gw, i, factors+[(gn,gw)], val*gv)
    rec(w, 0, [], mpf(1))
    return out
def pslq_identify(target, w, maxcoeff=10**6):
    mons = monomials(w)
    basis = [v for _,v in mons]
    rel = mpmath.pslq(basis + [target], tol=mpf(10)**(-70), maxcoeff=maxcoeff, maxsteps=10000)
    if rel is None: return None, mons
    # relation: sum rel[i]*basis[i] + rel[-1]*target = 0
    if rel[-1]==0: return None, mons
    return rel, mons
def show(rel, mons, target):
    terms = []
    for r,(f,v) in zip(rel, mons):
        if r:
            name = '*'.join(g for g,_ in f) if f else '1'
            terms.append(f"({r})*{name}")
    return ' + '.join(terms) + f" + ({rel[-1]})*T = 0"

ck = pickle.load(open('ckpt_integrals90.pkl','rb'))
I  = {eval(k): mpf(v) for k,v in ck['I'].items()}
results = {}
for key in [(0,2),(1,2),(0,3),(2,2),(1,3),(0,4),(3,2),(2,3),(1,4),(0,5)]:
    a,b = key; w = a+b+1
    rel, mons = pslq_identify(I[key], w)
    txt = show(rel,mons,I[key]) if rel else "NONE"
    results[str(key)] = txt
    print(f"I{key} w{w}: {txt}")

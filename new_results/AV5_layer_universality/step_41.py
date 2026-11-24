import mpmath, pickle, json
from mpmath import mp, mpf, zeta, log, polylog, nstr
mp.dps = 120
def make_gens():
    return [(1,'L',log(2)),(2,'Z2',zeta(2)),(3,'Z3',zeta(3)),(4,'Li4',polylog(4,mpf('0.5'))),
            (5,'Z5',zeta(5)),(5,'Li5',polylog(5,mpf('0.5'))),(6,'Li6',polylog(6,mpf('0.5'))),
            (7,'Z7',zeta(7)),(7,'Li7',polylog(7,mpf('0.5'))),(8,'Li8',polylog(8,mpf('0.5'))),
            (9,'Z9',zeta(9)),(9,'Li9',polylog(9,mpf('0.5')))]
gen = make_gens()
def monomials(w):
    out = []
    def rec(rem, start, factors, val):
        if rem==0: out.append((factors, val)); return
        for i in range(start, len(gen)):
            gw, gn, gv = gen[i]
            if gw <= rem: rec(rem-gw, i, factors+[gn], val*gv)
    rec(w, 0, [], mpf(1))
    return out
def identify(target, w, dps_for_bound=115):
    mons = monomials(w)
    N = len(mons)+1
    rel = mpmath.pslq([v for _,v in mons]+[target], maxcoeff=10**12, maxsteps=20000)
    if rel is None or rel[-1]==0: return None, mons, N
    return rel, mons, N
def fmt(rel, mons):
    terms=[]
    for r,(f,v) in zip(rel,mons):
        if r:
            terms.append(f"({r})" + ("*" + "*".join(f) if f else ""))
    return " + ".join(terms) + f" + ({rel[-1]})*T"
ck = pickle.load(open('ckpt_integrals135.pkl','rb'))
I  = {eval(k): mpf(v) for k,v in ck['I'].items()}
targets = [(0,2),(1,2),(0,3),(2,2),(1,3),(0,4),(3,2),(2,3),(1,4),(0,5),
           (4,2),(3,3),(2,4),(1,5),(0,6)]
out = {}
for key in targets:
    a,b=key; w=a+b+1
    rel,mons,N = identify(I[key], w)
    if rel:
        mc = max(abs(x) for x in rel)
        bound = mpf(10)**(mpf(115)/(N-1))
        status = "GOOD" if mc < bound/10**6 else ("MARGINAL" if mc < bound/100 else "SPURIOUS?")
        out[str(key)] = (fmt(rel,mons), str(mc), status, w)
        print(f"I{key} w{w} N{N} maxcoeff {nstr(mc,4)} pigeon {nstr(bound,4)} [{status}]")
        print("   ", fmt(rel,mons)[:150])
    else:
        out[str(key)] = ("NONE", "", "OPEN", w)
        print(f"I{key} w{w} N{N}: NONE")
json.dump(out, open('pslq_I_120.json','w'), indent=1)

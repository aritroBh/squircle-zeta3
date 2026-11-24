import mpmath, pickle
from mpmath import mp, mpf, zeta, log, polylog, nstr
mp.dps = 120
gen = [(1,'L',log(2)),(2,'Z2',zeta(2)),(3,'Z3',zeta(3)),(4,'Li4',polylog(4,mpf('0.5'))),
       (5,'Z5',zeta(5)),(5,'Li5',polylog(5,mpf('0.5'))),(6,'Li6',polylog(6,mpf('0.5'))),
       (7,'Z7',zeta(7)),(7,'Li7',polylog(7,mpf('0.5'))),(8,'Li8',polylog(8,mpf('0.5'))),
       (9,'Z9',zeta(9)),(9,'Li9',polylog(9,mpf('0.5')))]
def monomials(w):
    out = []
    def rec(rem, start, factors, val):
        if rem==0: out.append((factors, val)); return
        for i in range(start, len(gen)):
            gw, gn, gv = gen[i]
            if gw <= rem: rec(rem-gw, i, factors+[gn], val*gv)
    rec(w, 0, [], mpf(1))
    return out
def fmt(rel, mons):
    terms=[]
    for r,(f,v) in zip(rel,mons):
        if r: terms.append(f"({r})"+("*"+"*".join(f) if f else ""))
    return " + ".join(terms) + f" + ({rel[-1]})*T"
F = {k: mpf(v) for k,v in pickle.load(open('ckpt_Fk_125.pkl','rb')).items()}
for k in range(3,9):
    w = k+1
    mons = monomials(w)
    rel = mpmath.pslq([v for _,v in mons]+[F[k]], maxcoeff=10**10, maxsteps=20000)
    if rel and rel[-1]!=0:
        mc = max(abs(x) for x in rel); N=len(mons)+1
        bound = mpf(10)**(mpf(115)/(N-1))
        status = "GOOD" if mc < bound/10**6 else ("MARGINAL" if mc < bound/100 else "SPURIOUS")
        print(f"F_{k} w{w} N{N} maxc {nstr(mc,4)} bound {nstr(bound,3)} [{status}]")
        print("   ", fmt(rel,mons)[:200])
    else:
        print(f"F_{k} w{w}: NONE")

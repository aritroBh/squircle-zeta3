import mpmath, pickle
from mpmath import mp, mpf, zeta, log, polylog, nstr
mp.dps = 190
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
mons8 = monomials(8)
print("w8 classical basis size:", len(mons8))
ckB = pickle.load(open('ckpt_integrals210.pkl','rb'))
I  = {eval(k): mpf(v) for k,v in ckB['I'].items()}
K2 = {eval(k): mpf(v) for k,v in ckB['K2'].items()}
targets = {f"I{ k}":v for k,v in I.items() if k in [(4,3),(2,5),(0,7),(5,2),(1,6),(3,4),(0,8)]}
targets["K2(6)"] = K2[6]
for name, T in targets.items():
    rel = mpmath.pslq([v for _,v in mons8]+[T], maxcoeff=10**10, maxsteps=50000)
    if rel and rel[-1]!=0:
        mc = max(abs(x) for x in rel)
        print(f"{name}: relation maxc {nstr(mc,4)} (pigeonhole ~1e7) -- {'candidate' if mc<10**5 else 'spurious?'}")
    else:
        print(f"{name}: NONE")

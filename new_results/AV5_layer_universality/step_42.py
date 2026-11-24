import mpmath, pickle, json
from mpmath import mp, mpf, zeta, log, polylog, nstr
mp.dps = 125
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
def parse(relstr):
    # parse "(a)*X*Y + (b)*Z + (c)*T" format back into (coeffs, names)
    import re
    parts = relstr.split(" + ")
    out=[]
    for p in parts:
        m = re.match(r"\((-?\d+)\)(\*?(.*))", p)
        c = int(m.group(1)); rest = m.group(3)
        out.append((c, rest))
    return out
val_of = {gn: gv for _,gn,gv in gen}
val_of[''] = mpf(1)
data = json.load(open('pslq_I_120.json'))
ck = pickle.load(open('ckpt_integrals135.pkl','rb'))
I  = {eval(k): mpf(v) for k,v in ck['I'].items()}
print(f"{'target':>10} {'weight':>6} {'residual':>12} {'maxc*1e-124':>14} verdict")
for key,(relstr, mc, status, w) in data.items():
    if relstr=="NONE": continue
    terms = parse(relstr)
    R = mpf(0); maxc = 0
    for c, expr in terms:
        maxc = max(maxc, abs(c))
        if expr.endswith("*T") or expr=="T":
            targ = I[eval(key)]
            R += c*targ
        else:
            v = mpf(1)
            for f in expr.split("*"):
                v *= val_of[f]
            R += c*v
    ok = abs(R) < maxc*mpf(10)**(-123)
    print(f"{key:>10} {w:>6} {nstr(abs(R),4):>12} {nstr(maxc*mpf(10)**(-123),4):>14} {'VERIFIED' if ok else 'REJECTED'}")

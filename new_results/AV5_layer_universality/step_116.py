import mpmath, json
from mpmath import mp, mpf, sqrt, pi, zeta, log, nstr, matrix, lu_solve, catalan
mp.dps = 50
epss = [mpf(s) for s in json.load(open('ckpt_oct4_epss.json'))]
Ps = [mpf(s) for s in json.load(open('ckpt_oct4_Ps.json'))]
Poct = 16*(sqrt(2)-1)
def extract(taus, Ps, P0, lead, nterms):
    taus, Ps = taus[:nterms], Ps[:nterms]
    Q = [(P-P0)/t**lead for P,t in zip(Ps,taus)]
    A = matrix(nterms,nterms); b = matrix(nterms,1)
    for i,t in enumerate(taus):
        for j in range(nterms): A[i,j] = t**j
        b[i] = Q[i]
    return lu_solve(A,b)
c = extract(epss, Ps, Poct, 1, 6)
for j,ci in enumerate(c):
    print(f"d{j+1} = {nstr(ci,40)}")
L2 = log(2)
cand = {'1':mpf(1),'pi':pi,'pi^2':pi**2,'L2':L2,'L2^2':L2**2,'z2':zeta(2),'z3':zeta(3),'C':catalan,
        'sq2':sqrt(2),'sq2pi':sqrt(2)*pi,'sq2L2':sqrt(2)*L2,'piL2':pi*L2}
for j,T in enumerate(c[:3]):
    basis = list(cand.values()); names=list(cand.keys())
    rel = mpmath.pslq(basis+[T], maxcoeff=10**7)
    if rel and rel[-1]!=0:
        mc=max(abs(x) for x in rel)
        txt = " + ".join(f"({r})*{nm}" for r,nm in zip(rel,names) if r) + f" + ({rel[-1]})*d{j+1}"
        print(f"d{j+1}: maxc {nstr(mc,3)}: {txt[:160]}")
    else:
        print(f"d{j+1}: NONE")
# pinning check at the edge normal theta=0 (edge (1,-(√2-1)) to (1,√2-1): normal 0 deg)
def r_oct(phi, p):
    cc, ss = cos(phi), sin(phi)
    return (cc**p + ss**p + mpf(2)**(-p/2)*((cc+ss)**p + (cc-ss)**p))**(-1/p)
for p in [mpf(100), mpf(1000)]:
    # support at theta=0 vs octagon h(0)=1: r_oct(0) = (1+0+2^{-p/2}(1+1))^{-1/p} = (1+2^{1-p/2})^{-1/p}
    r0 = r_oct(0, p)
    print(f"p={p}: h(0) = r(0) = {nstr(r0,30)}, 1-r0 = {nstr(1-r0,8)}")

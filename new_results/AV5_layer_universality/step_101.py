import mpmath, json
from mpmath import mp, mpf, sqrt, pi, zeta, log, nstr, matrix, lu_solve, catalan
mp.dps = 70
epss = [mpf(s) for s in json.load(open('ckpt_hex3_epss.json'))]
Ps = [mpf(s) for s in json.load(open('ckpt_hex3_Ps.json'))]
Phex = 2*sqrt(2)+4
def extract(taus, Ps, P0, lead, nterms):
    taus, Ps = taus[:nterms], Ps[:nterms]
    n = nterms
    Q = [(P-P0)/t**lead for P,t in zip(Ps,taus)]
    A = matrix(n,n); b = matrix(n,1)
    for i,t in enumerate(taus):
        for j in range(n): A[i,j] = t**j
        b[i] = Q[i]
    return lu_solve(A,b)
# P = Phex + d1 eps + d2 eps^2 + ... (d1<0)
c = extract(epss, Ps, Phex, 1, 7)
for j,ci in enumerate(c):
    print(f"d{j+1} = {nstr(ci,40)}")
# PSLQ d1, d2
L2 = log(2)
cand = {'1':mpf(1),'pi':pi,'pi^2':pi**2,'pi^3':pi**3,'L2':L2,'L2^2':L2**2,'z2':zeta(2),'z3':zeta(3),
        'C':catalan,'piL2':pi*L2,'sq2':sqrt(2),'sq2pi':sqrt(2)*pi,'sq2L2':sqrt(2)*L2,'z2sq2':zeta(2)*sqrt(2)}
for j,T in enumerate(c[:3]):
    basis = list(cand.values()); names=list(cand.keys())
    rel = mpmath.pslq(basis+[T], maxcoeff=10**7)
    if rel and rel[-1]!=0:
        mc=max(abs(x) for x in rel)
        txt = " + ".join(f"({r})*{nm}" for r,nm in zip(rel,names) if r) + f" + ({rel[-1]})*d{j+1}"
        print(f"d{j+1}: maxc {nstr(mc,3)}: {txt[:170]}")
    else:
        print(f"d{j+1}: NONE")

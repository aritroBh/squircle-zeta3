from mpmath import mp, mpf, pi, log, zeta, pslq, catalan, nstr
mp.dps = 55
import json
res = json.load(open('scratchpad/bnder60.json'))
ln2 = log(2); G = catalan
ns=[1,3,5,7,9,11,13]
def show(n,k,basis,names,tol=mpf('1e-42'),mc=10**7):
    v = mpf(res[f"{n},{k}"])
    r = pslq([v]+basis, tol=tol, maxcoeff=mc)
    if r is None:
        print(f"n={n} k={k}: NO RELATION")
        return None
    # residual
    resid = r[0]*v + sum(a*b for a,b in zip(r[1:],basis))
    terms = [f"({r[0]})*b"]
    for a,nm in zip(r[1:],names):
        if a: terms.append(f"{a:+d}*{nm}")
    print(f"n={n} k={k}: {' '.join(terms)}  resid={nstr(resid,4)}")
    return r

print("=== k=2 basis {1,ln2,ln^2 2,pi^2,pi^2 ln2} ===")
b2 = [mpf(1),ln2,ln2**2,pi**2,pi**2*ln2]
nm2 = ["1","ln2","ln2^2","pi^2","pi^2 ln2"]
for n in ns: show(n,2,b2,nm2)

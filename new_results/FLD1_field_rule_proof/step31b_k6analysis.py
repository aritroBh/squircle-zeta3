"""step31b_k6analysis.py -- k=6 properness + weight budget from saved mu_zsum_K6.pkl."""
import pickle, sympy as sp
n=sp.Symbol('n')
Z=pickle.load(open('mu_zsum_K6.pkl','rb'))   # {zkey: [(w, cstr), ...]}
def wt(zkey): return sum((i+2)*e for i,e in enumerate(zkey))
def analyze(zkey, terms):
    maxdeg=-999; poles=set(); worst=0
    for w,cs in terms:
        e=sp.together(sp.sympify(cs)); num,den=sp.fraction(e)
        dd=sp.Poly(sp.expand(den),n); dn=sp.Poly(sp.expand(num),n)
        maxdeg=max(maxdeg, dn.degree()-dd.degree())
        wl = sum(abs(x) for x in w)
        if dd.degree()>0:
            rts=sp.roots(sp.expand(den),n)
            for rt,mult in rts.items(): poles.add(sp.nsimplify(rt))
            worst=max(worst, wl+max(rts.values()))
        else:
            worst=max(worst, wl)
    return maxdeg, poles, worst
allpoles=set(); allok=True
for zkey,terms in sorted(Z.items()):
    md,poles,worst=analyze(zkey,terms)
    allpoles|=poles
    bound=6+1-wt(zkey)
    okb = worst<=bound and md<=-2
    allok = allok and okb
    print(f"  z={zkey} (wt {wt(zkey)}): deg={md}, max(|w|+po)={worst} (bound {bound}) ok={okb}")
print("k=6: all poles:", sorted(allpoles, key=sp.N))
print("k=6 proper (deg<=-2) AND weight budget (|w|+po <= k+1-wt(z)):", allok)

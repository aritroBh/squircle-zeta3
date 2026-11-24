import sympy as sp
r3 = sp.sqrt(3)
Xs = {1: 1-2*r3/3, 2: 1-r3, 5: 1-r3/3, 6: 1-r3/2}
IDX=[1,2,5,6]
gam = {}; dl = {}
for a in IDX:
    for k in IDX:
        if a==k: continue
        gam[(a,k)] = sp.simplify(Xs[k]/(Xs[k]-Xs[a]))
        dl[(a,k)]  = sp.simplify(gam[(a,k)]*(1-Xs[a]))
extra = set()
for (a,k) in gam:
    extra.add(sp.simplify((gam[(a,k)]-1)/gam[(a,k)]))
    extra.add(sp.simplify((dl[(a,k)]-1)/dl[(a,k)]))
gd = set(gam.values())|set(dl.values())
print("|gamma/delta set| =", len(gd), "|extra set| =", len(extra))
for e in sorted(extra, key=lambda e: float(sp.N(e,10))):
    print("   ", e, "~", sp.N(e,7), "in gd?" , e in gd)

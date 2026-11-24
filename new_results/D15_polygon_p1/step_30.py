import mpmath as mp
mp.mp.dps = 170
P3_3 = mp.mpf('-5.145455017135615837171681131123975071115553566419741720927607420040126819690946561155264766767817630762120256954983911600508189884511766926')
sq3 = mp.sqrt(3)
l2 = mp.log(2); l3 = mp.log(3); pi = mp.pi
Cl = mp.clsin(2, pi/3); Li213 = mp.polylog(2, mp.mpf(1)/3)
Li313 = mp.polylog(3, mp.mpf(1)/3); Li312 = mp.polylog(3, mp.mpf(1)/2)
z3 = mp.zeta(3)
T = {'1':mp.mpf(1),'l2':l2,'l3':l3,'pi':pi,'pi2':pi**2,'l2^2':l2**2,'l3^2':l3**2,'l2l3':l2*l3,
     'Cl':Cl,'Li213':Li213,'pil2':pi*l2,'pil3':pi*l3,
     'z3':z3,'pi3':pi**3,'pi2l2':pi**2*l2,'pi2l3':pi**2*l3,'pil2^2':pi*l2**2,'pil3^2':pi*l3**2,
     'pil2l3':pi*l2*l3,'l2^3':l2**3,'l3^3':l3**3,'l2^2l3':l2**2*l3,'l2l3^2':l2*l3**2,
     'Cll2':Cl*l2,'Cll3':Cl*l3,'piCl':pi*Cl,'Li313':Li313,'Li312':Li312,
     'Li213l2':Li213*l2,'Li213l3':Li213*l3,'piLi213':pi*Li213}
names = list(T.keys()); tvals = [T[n] for n in names]
print("basis size:", len(names)+1)
target = P3_3/sq3
rel = mp.pslq([target]+tvals, tol=mp.mpf('1e-110'), maxcoeff=10**6, maxsteps=50000)
if rel is None:
    print("NO relation at maxcoeff 1e6")
else:
    resid = target*rel[0]+sum(r*v for r,v in zip(rel[1:],tvals))
    print("relation, residual:", mp.nstr(resid,5), "maxcoeff:", max(abs(r) for r in rel))
    for n,r in zip(['P3/sq3']+names, rel):
        if r: print(f"  {r} * {n}")

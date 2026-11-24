import mpmath as mp
mp.mp.dps = 100
P2_3 = mp.mpf('-1.00109443973840992962135105887535450171338492998892655629450673110489944312676081385506309')
sq3 = mp.sqrt(3)
l2 = mp.log(2); l3 = mp.log(3)
Cl = mp.clsin(2, mp.pi/3)
Li2_13 = mp.polylog(2, mp.mpf(1)/3)
T = {'1':mp.mpf(1), 'pi':mp.pi, 'l2':l2, 'l3':l3, 'pi2':mp.pi**2, 'pi*l2':mp.pi*l2, 'pi*l3':mp.pi*l3,
     'l2^2':l2**2, 'l3^2':l3**2, 'l2*l3':l2*l3, 'Li2(1/3)':Li2_13, 'Cl':Cl}
names = list(T.keys()); tvals = [T[n] for n in names]
basis = []; bn = []
for f,fn in [(mp.mpf(1),'1'),(sq3,'sq3')]:
    for n,v in zip(names,tvals):
        basis.append(f*v); bn.append(f'{fn}*{n}')
rel = mp.pslq([P2_3]+basis, tol=mp.mpf('1e-70'), maxcoeff=10**10, maxsteps=20000)
if rel is None:
    print("no relation in 12x2 basis; trying smaller cores")
else:
    resid = P2_3*rel[0]+sum(r*b for r,b in zip(rel[1:],basis))
    print("relation found, residual:", mp.nstr(resid,5), "maxcoeff:", max(abs(r) for r in rel))
    for n,r in zip(['P2']+bn, rel):
        if r: print(f"  {r} * {n}")
# also try core bases
for core in [['1','pi2'], ['1','pi2','l2^2','l3^2','l2*l3'], ['pi2','l2^2','l3^2','l2*l3','Li2(1/3)','Cl']]:
    basis2 = []; bn2 = []
    for f,fn in [(mp.mpf(1),'1'),(sq3,'sq3')]:
        for n in core:
            basis2.append(f*T[n]); bn2.append(f'{fn}*{n}')
    rel = mp.pslq([P2_3]+basis2, tol=mp.mpf('1e-75'), maxcoeff=10**10, maxsteps=20000)
    print(core, "->", rel)

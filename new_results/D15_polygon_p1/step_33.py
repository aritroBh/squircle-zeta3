import mpmath as mp
mp.mp.dps = 170
P2_4 = mp.mpf('-1.4562051103812046446537306966795936456777152413012283082812637083293482631402138804126453990677125615381393921189383069887917998883431152326')
sq2 = mp.sqrt(2); c8 = mp.cos(mp.pi/8); s8 = mp.sin(mp.pi/8)
l2 = mp.log(2); L = mp.log(1+sq2); pi = mp.pi
T = {'1':mp.mpf(1), 'l2':l2, 'L':L, 'pi2':pi**2, 'l2^2':l2**2, 'L^2':L**2, 'l2L':l2*L}
names=list(T.keys()); tvals=[T[n] for n in names]
for label, target in [('P2/c8', P2_4/c8), ('P2/(sq2c8)', P2_4/(sq2*c8)), ('P2/s8', P2_4/s8), ('P2/(sq2s8)', P2_4/(sq2*s8)), ('P2/P0', P2_4/mp.mpf('2.53620268449771540274021154137584885279564001427821445195994336040708703329476758785299783'))]:
    rel = mp.pslq([target]+tvals, tol=mp.mpf('1e-120'), maxcoeff=10**8, maxsteps=50000)
    if rel is None:
        print(label, ": none")
        continue
    resid = target*rel[0]+sum(r*v for r,v in zip(rel[1:],tvals))
    mc = max(abs(r) for r in rel)
    print(label, "residual:", mp.nstr(resid,4), "maxcoeff:", mc, " rel:", rel)

import mpmath as mp
mp.mp.dps = 170
P2_4 = mp.mpf('-1.4562051103812046446537306966795936456777152413012283082812637083293482631402138804126453990677125615381393921189383069887917998883431152326')
sq2 = mp.sqrt(2); c8 = mp.cos(mp.pi/8)
l2 = mp.log(2); L = mp.log(1+sq2); pi = mp.pi
field4 = {'1':mp.mpf(1), 'sq2':sq2, 'c8':c8, 'sq2c8':sq2*c8}
T = {'1':mp.mpf(1), 'l2':l2, 'L':L, 'pi2':pi**2, 'l2^2':l2**2, 'L^2':L**2, 'l2L':l2*L}
basis = []; bn = []
for fn,f in field4.items():
    for tn,t in T.items():
        basis.append(f*t); bn.append(f'{fn}*{tn}')
rel = mp.pslq([P2_4]+basis, tol=mp.mpf('1e-110'), maxcoeff=10**6, maxsteps=50000)
if rel is None:
    print("no relation; extend basis with Clausen/Li2")
else:
    resid = P2_4*rel[0]+sum(r*v for r,v in zip(rel[1:],basis))
    print("relation, residual:", mp.nstr(resid,5), "maxcoeff:", max(abs(r) for r in rel), " rel[0]:", rel[0])
    for n,r in zip(['P2']+bn, rel):
        if r: print(f"  {r} * {n}")

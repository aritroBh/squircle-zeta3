import mpmath as mp
mp.mp.dps = 100
P1_4 = mp.mpf('2.75023713004966371870863005296385194252854498484351236290611395868134254516571310482650575')
P1_3 = mp.mpf('2.40113226770588735429020690407746529554651722044431670825375284904281178764684126665661951')
# m=4 PSLQ: field Q(cos(pi/8)) x {1, log2, log(1+sqrt2)}
c8 = mp.cos(mp.pi/8); sq2 = mp.sqrt(2)
field4 = [mp.mpf(1), sq2, c8, sq2*c8]
trans4 = [mp.mpf(1), mp.log(2), mp.log(1+sq2)]
basis = [f*t for f in field4 for t in trans4]
rel = mp.pslq([P1_4] + basis, tol=mp.mpf('1e-80'), maxcoeff=10**12, maxsteps=10000)
print("m=4 P' relation:", rel)
if rel:
    resid = P1_4*rel[0] + sum(r*b for r,b in zip(rel[1:], basis))
    print("residual:", mp.nstr(resid, 5), " maxcoeff:", max(abs(r) for r in rel))
# m=3 P' recheck with pslq too
sq3 = mp.sqrt(3)
rel3 = mp.pslq([P1_3, sq3*mp.log(2)], tol=mp.mpf('1e-85'), maxcoeff=10**12, maxsteps=10000)
print("m=3 P' vs sq3*log2:", rel3, "residual:", mp.nstr(P1_3*rel3[0]+sq3*mp.log(2)*rel3[1],5))
# m=4: also try ratio P'/P(1)
P0_4 = mp.mpf('2.53620268449771540274021154137584885279564001427821445195994336040708703329476758785299783')
print("P'/P(1) m=4:", mp.nstr(P1_4/P0_4, 50))

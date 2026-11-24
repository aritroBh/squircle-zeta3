import mpmath as mp
mp.mp.dps = 100
P1_3 = mp.mpf('2.40113226770588735429020690407746529554651722044431670825375284904281178764684126665661951')
sq3 = mp.sqrt(3)
print("P'/sqrt3 =", mp.nstr(P1_3/sq3, 60))
print("log4     =", mp.nstr(mp.log(4), 60))
print("diff     =", mp.nstr(P1_3/sq3 - mp.log(4), 10))
print("P'/(2sqrt3) =", mp.nstr(P1_3/(2*sq3), 60))
print("log2        =", mp.nstr(mp.log(2), 60))
print("diff =", mp.nstr(P1_3/(2*sq3) - mp.log(2), 10))
# quick PSLQ attempts for P'(1) m=3
Cl = mp.clsin(2, mp.pi/3)
print("Cl2(pi/3) =", mp.nstr(Cl, 60))
for basis_name, basis in [
    ("A: 1,sq3,log2,log3", [mp.mpf(1), sq3, mp.log(2), mp.log(3)]),
    ("B: +Cl,pi", [mp.mpf(1), sq3, mp.log(2), mp.log(3), Cl, mp.pi]),
    ("C: +sq3*log2,sq3*log3,sq3*Cl,sq3*pi", [mp.mpf(1), sq3, mp.log(2), mp.log(3), Cl, mp.pi,
                                              sq3*mp.log(2), sq3*mp.log(3), sq3*Cl, sq3*mp.pi]),
]:
    rel = mp.pslq([P1_3] + basis, tol=mp.mpf('1e-85'), maxcoeff=mp.mpf('1e12'), maxsteps=10000)
    print(basis_name, "->", rel)

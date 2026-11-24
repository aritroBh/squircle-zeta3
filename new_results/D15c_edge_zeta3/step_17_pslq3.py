import mpmath as mp
mp.mp.dps=150
import json
vals=json.load(open('vals_m2_8_200.json'))
z3=mp.zeta(3);pi=mp.pi
def Cm(m):
    h=pi/(2*m); s=mp.sin(h); return 48*m*s**4*mp.cos(h)*(3*s**2-1)
# m=2 known full: 6sqrt2 z3 - 2sqrt2 pi^2
p2=mp.mpf(vals['2'][3])
print('m=2: P3 - (6sqrt2 z3 - 2sqrt2 pi^2) resid:', mp.nstr(p2-(6*mp.sqrt(2)*z3-2*mp.sqrt(2)*pi**2),3))
# m=3: PSLQ against z3 + z3-free basis
p3=mp.mpf(vals['3'][3])
sq3=mp.sqrt(3)
# z3-free weight<=3 basis (algebraic x consts), for hexagon field Q(sqrt3)
bfree=[sq3*pi**3, pi**3, sq3*pi**2, pi**2, sq3*pi**2*mp.log(3), pi**2*mp.log(3), sq3*pi*mp.log(3)**2, pi*mp.log(3)**2, sq3*mp.log(3)**3, mp.log(3)**3, sq3*pi**2*mp.log(2), pi**2*mp.log(2)]
bfn=['sq3pi^3','pi^3','sq3pi^2','pi^2','sq3pi^2log3','pi^2log3','sq3pi log3^2','pi log3^2','sq3 log3^3','log3^3','sq3pi^2log2','pi^2log2']
rel=mp.pslq(mp.matrix([p3, z3]+bfree), tol=mp.mpf('1e-120'), maxcoeff=10**6)
print('m=3 pslq:', rel)
if rel is not None and rel[0]!=0:
    print('   z3 coeff of P3 =', mp.nstr(-mp.mpf(rel[1])/rel[0],10), '  C_3=', mp.nstr(Cm(3),10))
    # reconstruct
    val = -mp.mpf(rel[1])/rel[0]*z3 - sum(mp.mpf(rel[k+2])/rel[0]*bfree[k] for k in range(len(bfree)))
    print('   resid:', mp.nstr(p3-val,3))


import mpmath as mp, json
mp.mp.dps = 110
d = json.load(open("m3_total.json"))
rem = mp.mpf(d["rem"])
sqrt3=mp.sqrt(3); pi=mp.pi; l2=mp.log(2); l3=mp.log(3); z3=mp.zeta(3); G=mp.catalan
# independent elementary basis (drop z2 = pi^2/6). weight<=5 hexagon constants.
basis = [mp.mpf(1), sqrt3, pi, pi*sqrt3, pi**2, pi**2*sqrt3, pi**3, pi**3*sqrt3,
         pi**4, pi**4*sqrt3, pi**5, pi**5*sqrt3,
         l2, l3, l2*sqrt3, l3*sqrt3, l2**2, l3**2, l2*l3, l3**2*sqrt3,
         z3, z3*sqrt3, z3*pi, z3*pi*sqrt3, z3*l2, z3*l3,
         G, G*sqrt3, G*pi, G*pi*sqrt3, l3**3, l2**3]
bn = ['1','r3','pi','pir3','pi2','pi2r3','pi3','pi3r3','pi4','pi4r3','pi5','pi5r3',
      'l2','l3','l2r3','l3r3','l2^2','l3^2','l2l3','l3^2r3',
      'z3','z3r3','z3pi','z3pir3','z3l2','z3l3','G','Gr3','Gpi','Gpir3','l3^3','l2^3']
rel = mp.pslq(mp.matrix([rem]+basis), tol=mp.mpf('1e-92'), maxcoeff=10**12, maxsteps=50000)
print("PSLQ:", rel)
if rel is not None and rel[0]!=0:
    rel=list(rel)
    print("rem =", {bn[k]: mp.nstr(-mp.mpf(rel[k+1])/rel[0],10) for k in range(len(basis)) if rel[k+1]!=0})
else:
    print("no relation with rem involved; trying to guess structure")
    # Try rem against pi^5*sqrt3 etc individual
    for nm,b in zip(bn,basis):
        print("  rem/%s = %s"%(nm, mp.nstr(rem/b,12)))

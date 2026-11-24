import mpmath as mp
mp.mp.dps = 50
for m in [2,3,4]:
    h = mp.pi/(2*m); s = mp.sin(h); c = mp.cos(h)
    C5 = 2880*m*s**4*c*(65*s**6 - 63*s**4 + 17*s**2 - 1)
    C23 = -960*m*c*s**4*(45*s**6 - 27*s**4 + 6*s**2 - 1)
    print(f"m={m}: formula C5 = {mp.nstr(C5, 25)}")
    print(f"      formula C23 = {mp.nstr(C23, 25)}")
print("D17 pert: m=3 C5 = 153.44887623305522272  (resid 2e-59 vs their closed)")
print("my m=3 formula:", mp.nstr(2880*3*mp.sin(mp.pi/6)**4*mp.cos(mp.pi/6)*(65*mp.sin(mp.pi/6)**6-63*mp.sin(mp.pi/6)**4+17*mp.sin(mp.pi/6)**2-1), 25))
# maybe their polygon is normalized differently (e.g. apothem 1, R=1)? check ratio:
val = 153.44887623305522272
mine = 2880*3*mp.sin(mp.pi/6)**4*mp.cos(mp.pi/6)*(65*mp.sin(mp.pi/6)**6-63*mp.sin(mp.pi/6)**4+17*mp.sin(mp.pi/6)**2-1)
print("ratio D17/mine:", mp.nstr(mp.mpf(val)/mine, 25))

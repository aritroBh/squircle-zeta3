import os, mpmath as mp
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 100
phi=(1+mp.sqrt(5))/2; lf=mp.log(phi); z2=mp.pi**2/6; z3=mp.zeta(3)
v3 = mp.polylog(3,-phi)
print("Li3(-phi):", repr(v3))
v2 = mp.polylog(3,-2)
print("Li3(-2):", repr(v2))
# derived inversion: F(x) = Li3(-x)+Li3(-1/x) = 2 Li3(-1/x) - z2 ln x - ln^3 x/6
for x in [phi, mp.mpf(2), mp.mpf(3)]:
    F = mp.polylog(3,-x) + mp.polylog(3,-1/x)
    cand = 2*mp.polylog(3,-1/x) - z2*mp.log(x) - mp.log(x)**3/6
    print(f"x={mp.nstr(x,8)}: F={mp.nstr(F,20)}  cand={mp.nstr(cand,20)}  diff={mp.nstr(abs(F-cand),6)}")
# alternative candidate with +z2 ln x? derive again but check numerically for pattern: F + z2 lnx + ln^3 x/6 vs 2 Li3(-1/x)
for x in [phi, mp.mpf(2)]:
    lhs = mp.polylog(3,-x) - mp.polylog(3,-1/x)
    print(f"x={mp.nstr(x,8)}: Li3(-x)-Li3(-1/x) = {mp.nstr(lhs,20)}; -z2 lnx - ln^3 x/6 = {mp.nstr(-z2*mp.log(x)-mp.log(x)**3/6,20)}")
    print("   sum check: ", mp.nstr(lhs + z2*mp.log(x) + mp.log(x)**3/6, 10))
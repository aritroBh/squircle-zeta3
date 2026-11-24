import mpmath as mp, time
mp.mp.dps = 120
s5 = mp.sqrt(5); phi = (1+s5)/2; lf = mp.log(phi); z2 = mp.pi**2/6; z3 = mp.zeta(3)
x1 = -1/phi; x4 = 1/phi**2
L2 = lambda x: mp.polylog(2,x); L3 = lambda x: mp.polylog(3,x)
def A(x): return mp.nsum(lambda m: x**m*mp.harmonic(m)/m**2, [1, mp.inf])
def Li21(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1); m=0
    while True:
        m+=1; xp*=x
        if m>=2: tot += xp*ps/m**2
        yp*=y; ps += yp/m
        if m>10 and abs(xp*ps/m**2) < mp.mpf('1e-125'): break
    return tot
Lam = Li21(x1,x1); swap = Li21(x4,-phi)
checks = []
checks.append(("golden Li2(1/phi^2) = 2z2/5 - lf^2", L2(x4) - (2*z2/5 - lf**2)))
checks.append(("golden Li2(1/phi)   = 3z2/5 - lf^2", L2(1/phi) - (3*z2/5 - lf**2)))
checks.append(("golden Li2(-1/phi)  = -2z2/5 + lf^2/2", L2(x1) - (-2*z2/5 + lf**2/2)))
checks.append(("AA(x1) = 4z3/5 - 2 lf z2/5 + 2 L3(x1)", A(x1) - (4*z3/5 - 2*lf*z2/5 + 2*L3(x1))))
checks.append(("AA(x4) = z3 - 3 lf z2/5 + 4 L3(x1) + 3 L3(1/phi)", A(x4) - (z3 - 3*lf*z2/5 + 4*L3(x1) + 3*L3(1/phi))))
checks.append(("Lam = 2z3/5 + 4 lf z2/5 - lf^3 + 3 L3(x1) + L3(1/phi)", Lam - (2*z3/5 + 4*lf*z2/5 - lf**3 + 3*L3(x1) + L3(1/phi))))
checks.append(("swap = 3z3/5 - 7 lf z2/5 + 2 lf^3 + L3(x1) + L3(1/phi)", swap - (3*z3/5 - 7*lf*z2/5 + 2*lf**3 + L3(x1) + L3(1/phi))))
checks.append(("L3(1/phi^2) = 4 L3(-1/phi) + 4 L3(1/phi) [Ra/duplication]", L3(x4) - (4*L3(x1)+4*L3(1/phi))))
checks.append(("Lewin: L3(1/phi^2) = 4(z3-lf z2)/5 + 2 lf^3/3", L3(x4) - (4*(z3-lf*z2)/5 + 2*lf**3/3)))
checks.append(("sum rel: L3(1/phi)+L3(-1/phi) = (z3-lf z2)/5 + lf^3/6", L3(1/phi)+L3(x1) - ((z3-lf*z2)/5 + lf**3/6)))
worst = mp.mpf(0)
for name, d in checks:
    print(f"{name:60s} |resid| = {mp.nstr(abs(d),3)}")
    worst = max(worst, abs(d))
print("worst:", mp.nstr(worst,3))
with open('<project>/research/08_theory/DEC2_decagon_exact_assembly/progress.log','a') as f:
    f.write(f"step_15: reduction statements validated at 120 dps (golden Li2 x3, AA x2 [DEC1-proved], Lam/swap [PRV1-proved], Ra, Lewin, sum rel): worst {mp.nstr(worst,2)}\n")

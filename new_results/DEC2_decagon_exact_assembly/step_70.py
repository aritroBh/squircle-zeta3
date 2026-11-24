import mpmath as mp
mp.mp.dps = 60
s5 = mp.sqrt(5); phi = (1+s5)/2; lf = mp.log(phi); z2 = mp.pi**2/6; z3 = mp.zeta(3)
Li3 = lambda x: mp.polylog(3,x)
# LHS of 3-term Landen at x = 1/phi^2
x = 1/phi**2
LHS = Li3(x) + Li3(1-x) + Li3(x/(x-1))
print("LHS (x=1/phi^2):", mp.nstr(LHS,40))
# candidate RHS: z3 + z2*log(1-x) - (1/2)log^2 x log(1-x) + c*log^3(1-x)
base = z3 + z2*mp.log(1-x) - mp.mpf('0.5')*mp.log(x)**2*mp.log(1-x)
print("base (without log^3 term):", mp.nstr(base,40))
print("residual/(log(1-x)^3):", mp.nstr((LHS-base)/mp.log(1-x)**3,20))
# also test at x=1/2
x2 = mp.mpf('0.5')
LHS2 = Li3(x2)+Li3(1-x2)+Li3(x2/(x2-1))
base2 = z3 + z2*mp.log(1-x2) - mp.mpf('0.5')*mp.log(x2)**2*mp.log(1-x2)
print("x=1/2: residual/log^3(1-x):", mp.nstr((LHS2-base2)/mp.log(1-x2)**3,20))
# and at x=0.3
x3 = mp.mpf('0.3')
LHS3 = Li3(x3)+Li3(1-x3)+Li3(x3/(x3-1))
base3 = z3 + z2*mp.log(1-x3) - mp.mpf('0.5')*mp.log(x3)**2*mp.log(1-x3)
print("x=0.3: residual/log^3(1-x):", mp.nstr((LHS3-base3)/mp.log(1-x3)**3,20))
# check DEC1 Lewin: Li3(1/phi^2) = 4(z3-lf*z2)/5 + 2 lf^3/3
print("Lewin check:", mp.nstr(Li3(1/phi**2) - (mp.mpf(4)/5*(z3-lf*z2) + mp.mpf(2)/3*lf**3),6))
# sum relation
print("sum rel check:", mp.nstr(Li3(1/phi)+Li3(-1/phi) - ((z3-lf*z2)/5 + lf**3/6),6))

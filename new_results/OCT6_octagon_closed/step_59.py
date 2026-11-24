import mpmath as mp
mp.mp.dps = 400
r2=mp.sqrt(2); k1=r2-1; k2=r2/2-1; x1=-k1; x2=-k2
l2=mp.log(2); lu=mp.log(1+r2); z2=mp.pi**2/6; z3=mp.zeta(3)
Li2=lambda z: mp.polylog(2,z); Li3=lambda z: mp.polylog(3,z)
def AA_of(xv):
    N = int(400/mp.log10(1/abs(xv)))+50
    s=mp.mpf(0); Hn=mp.mpf(0); xp=mp.mpf(1)
    for nn in range(1,N+1):
        Hn+=1/mp.mpf(nn); xp*=xv; s+=xp*Hn/nn**2
    return s
def AA_cont(xv):
    return mp.quad(lambda t: (mp.polylog(2,t)+mp.log(1-t)**2/2)/t, [0, xv])
# difference identity: A(x)-A(1-x) = 2Li3(x)-2Li3(1-x)+log(1-x)Li2(1-x)-logx Li2(x)-(1/2)log^2 x log(1-x)+(1/2)log^2(1-x) log x
def DIF(x):
    return (2*Li3(x)-2*Li3(1-x)+mp.log(1-x)*Li2(1-x)-mp.log(x)*Li2(x)
            -mp.log(x)**2*mp.log(1-x)/2+mp.log(1-x)**2*mp.log(x)/2)
for xx in [mp.mpf('0.333333'), x2, mp.mpf('0.117'), r2/2, mp.mpf('0.91')]:
    print("DIFF residual at x=%s:"%mp.nstr(xx,8), mp.nstr(abs(AA_of(xx)-AA_of(1-xx)-DIF(xx)),3))
# L2 and L1 RHS
def R_L2(x):
    return z2*mp.log(1-x) + 2*mp.log(1-x)*Li2(1-x) + mp.log(1-x)**2*mp.log(x)/2 - 3*Li3(1-x) + 3*z3
def R_L1(y):
    return -z2*mp.log(y) - mp.log(y)**3/3 + Li3(-y) + 2*Li3(-1/y) - mp.log(y)*Li2(-y) + z3
# solved A(x1): 2A(x1) = R_L2(x2) - DIF(x2) - R_L2(r2h) + R_L1(k1)
Ax1_sol = (R_L2(x2) - DIF(x2) - R_L2(r2/2) + R_L1(k1))/2
print("\nA(x1) analytic-solution =", mp.nstr(Ax1_sol,25))
print("A(x1) direct           =", mp.nstr(AA_of(x1),25))
print("|diff| =", mp.nstr(abs(Ax1_sol-AA_of(x1)),3))
# A(x2) too: A(x2) = R_L2(x2) - A(x1)
print("A(x2) via L2  =", mp.nstr(R_L2(x2)-Ax1_sol,25), " direct:", mp.nstr(AA_of(x2),25),
      " diff:", mp.nstr(abs(R_L2(x2)-Ax1_sol-AA_of(x2)),3))

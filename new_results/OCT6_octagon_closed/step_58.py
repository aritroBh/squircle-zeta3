import mpmath as mp, pickle
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
def R4t(x):
    return (2*Li3(1-x)-Li3(x)+mp.log(x)*Li2(x)-mp.log(x)*mp.log(1-x)**2/2+mp.log(x)**2*mp.log(1-x)/2
            + z2*mp.log(1-x) + 5*z3/2 - l2**3/6)
for xx in [mp.mpf('0.333333'), x2, mp.mpf('0.117'), r2/2]:
    lhs = AA_of(xx)+3*AA_of(1-xx)+AA_of(xx/(xx-1))
    print("4-term residual at x=%s:"%mp.nstr(xx,8), mp.nstr(abs(lhs-R4t(xx)),3))
def RL1(y):
    return -z2*mp.log(y) - mp.log(y)**3/3 + Li3(-y) + 2*Li3(-1/y) - mp.log(y)*Li2(-y) + z3
Am_r2 = AA_cont(-(1+r2))  # one continuation value, reuse
import time; t0=time.time()
Am_r2 = AA_cont(-(1+r2)); print("AA_cont(-(1+r2)) %.1fs"%(time.time()-t0), mp.nstr(Am_r2,25))
for yy in [mp.mpf('0.5'), mp.mpf('0.9'), mp.mpf('1.7'), k1]:
    lhs = AA_of(-yy) + (AA_cont(-1/yy) if abs(1/yy)>1 else AA_of(-1/yy))
    print("L1 residual at y=%s:"%mp.nstr(yy,8), mp.nstr(abs(lhs-RL1(yy)),3))
R_L2_x2 = (z2*mp.log(1-x2) + 2*mp.log(1-x2)*Li2(1-x2) + mp.log(1-x2)**2*mp.log(x2)/2 - 3*Li3(1-x2) + 3*z3)
Ax1_sol = (RL1(k1) - R4t(r2/2) + 3*R_L2_x2 + (R4t(x2)-R_L2_x2)/3)/4
print("A(x1) solved =", mp.nstr(Ax1_sol,25))
print("A(x1) direct =", mp.nstr(AA_of(x1),25))
print("|diff| =", mp.nstr(abs(Ax1_sol-AA_of(x1)),3))

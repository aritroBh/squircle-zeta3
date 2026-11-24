import mpmath as mp, pickle, sympy as sp
mp.mp.dps = 745
r2=mp.sqrt(2); k1=r2-1; x1=-k1; x2=-k2 if (k2:=r2/2-1) else None
l2=mp.log(2); lu=mp.log(1+r2); z2=mp.pi**2/6; z3=mp.zeta(3)
def AA_of(xv):
    N = int(750/mp.log10(1/min(abs(xv),mp.mpf('0.99'))))+50
    s=mp.mpf(0); Hn=mp.mpf(0); xp=mp.mpf(1)
    for nn in range(1,N+1):
        Hn+=1/mp.mpf(nn); xp*=xv; s+=xp*Hn/nn**2
    return s
Li2=lambda z: mp.polylog(2,z); Li3=lambda z: mp.polylog(3,z)
# Lemma L2 (Landen): A(x)+A(x/(x-1)) = z2 log(1-x) + 2 log(1-x) Li2(1-x) + (1/2) log^2(1-x) log x - 3 Li3(1-x) + 3 z3
def L2res(x):
    w = x/(x-1)
    return AA_of(x)+AA_of(w) - (z2*mp.log(1-x) + 2*mp.log(1-x)*Li2(1-x) + mp.log(1-x)**2*mp.log(x)/2 - 3*Li3(1-x) + 3*z3)
for xx in [mp.mpf('0.333333'), x2, mp.mpf('0.117')]:
    print("L2 residual at x=%s:"%mp.nstr(xx,6), mp.nstr(abs(L2res(xx)),3))
# Lemma L1: 3A(-x)+A(-1/x) = z2 log x + 3 Li3(-x) - log x Li2(-x) - z3/4   (x>0; analytic continuation for A<-1)
def L1res(x):
    return 3*AA_of(-x)+AA_of(-1/x) - (z2*mp.log(x) + 3*Li3(-x) - mp.log(x)*Li2(-x) - z3/4)
for xx in [mp.mpf('0.5'), mp.mpf('0.9'), mp.mpf('1.0')]:
    print("L1 residual at x=%s:"%mp.nstr(xx,6), mp.nstr(abs(L1res(xx)),3))
# A(x1) via elimination using L1 at k1 and at 1+sqrt2 (uses mpmath analytic continuation of A below -1)
def AA_cont(xv):
    # A via integral representation: A(x) = int_0^x [Li2(t)+1/2 log^2(1-t)]/t dt, mpmath quad path on negative axis
    f=lambda t: (mp.polylog(2,t)+mp.log(1-t)**2/2)/t
    return mp.quad(f, [0, xv])
Rk = z2*mp.log(k1) + 3*Li3(-k1) - mp.log(k1)*Li2(-k1) - z3/4
Rb = z2*mp.log(1+r2) + 3*Li3(-(1+r2)) - mp.log(1+r2)*Li2(-(1+r2)) - z3/4
Ax1_elim = (3*Rk - Rb)/8
Ax1_direct = AA_of(x1)
print("A(x1) direct   =", mp.nstr(Ax1_direct,25))
print("A(x1) via L1   =", mp.nstr(Ax1_elim,25))
print("|diff| =", mp.nstr(abs(Ax1_elim-Ax1_direct),3))
# check AA_cont consistency: A(-(1+r2)) continuation: L1 at k1: A(-(1+r2)) = R(k1)-3A(x1)
print("L1 at k1 residual (with continuation):", mp.nstr(abs(3*Ax1_direct + AA_cont(-(1+r2)) - Rk),3))

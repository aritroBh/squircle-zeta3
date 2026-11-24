# step_48: is L3(1/phi)-L3(-1/phi) = q*L(chi5,3) with rational q? serious pslq
import mpmath as mp
mp.mp.dps = 200
def Li3(x): return mp.re(mp.polylog(3, mp.mpf(x)))
def Li2(x): return mp.re(mp.polylog(2, mp.mpf(x)))
def L3(x):
    x = mp.mpf(x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(1-x)
def chi5(n):
    r = n % 5
    return 0 if r==0 else (1 if r in (1,4) else -1)
L5 = mp.nsum(lambda n: chi5(int(n))/mp.mpf(n)**3, [1, mp.inf])
z3 = mp.zeta(3)
phi = (1+mp.sqrt(5))/2
X1, Xm1 = L3(1/phi), L3(-1/phi)
D = X1 - Xm1
rel = mp.pslq([D, L5, z3], tol=mp.mpf('1e-170'), maxcoeff=10**30, maxsteps=10000)
print("pslq [D, L5, z3]:", rel)
# also test the odd part of {1/phi^2}: should be zero (Lewin)
Y1 = L3(1/phi**2); Ym1 = L3(phi**2) if phi**2 < 1 else L3(1/phi**2)
print("L3(1/phi^2) =", Y1, " vs (4/5)zeta3:", mp.mpf(4)/5*z3, " diff:", Y1-mp.mpf(4)/5*z3)
# test golden Li3 sum relation from DEC1: Li3(-1/phi)+Li3(1/phi) = (zeta3 - l_phi*pi^2/6)/5 + l_phi^3/6
lphi = mp.log(phi)
lhs = Li3(-1/phi)+Li3(1/phi)
rhs = (z3 - lphi*mp.pi**2/6)/5 + lphi**3/6
print("DEC1 sum relation diff:", lhs-rhs)
# individual Li3(1/phi), Li3(-1/phi) values
print("Li3(1/phi) =", Li3(1/phi))
print("Li3(-1/phi) =", Li3(-1/phi))
print("X1-Xm1 =", D)
print("D/L5 =", D/L5)

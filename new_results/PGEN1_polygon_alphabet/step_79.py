# step_51: correct L-values via Hurwitz zeta; redo golden beta and eq:chi3 control
import mpmath as mp
mp.mp.dps = 150
def Li3(x): return mp.re(mp.polylog(3, mp.mpf(x)))
def Li2(x): return mp.re(mp.polylog(2, mp.mpf(x)))
def L3(x):
    x = mp.mpf(x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(1-x)
z3 = mp.zeta(3)
# L(chi5,3) and L(chi8,3) via Hurwitz zeta
L5 = (mp.zeta(3, mp.mpf(1)/5) - mp.zeta(3, mp.mpf(2)/5) - mp.zeta(3, mp.mpf(3)/5) + mp.zeta(3, mp.mpf(4)/5))/125
L8 = (mp.zeta(3, mp.mpf(1)/8) - mp.zeta(3, mp.mpf(3)/8) - mp.zeta(3, mp.mpf(5)/8) + mp.zeta(3, mp.mpf(7)/8))/512
print("L(chi5,3) =", mp.nstr(L5, 60))
print("L(chi8,3) =", mp.nstr(L8, 60), " paper 0.95838045456...")
# direct series cross-check
def chi5(n):
    r = n % 5
    return 0 if r==0 else (1 if r in (1,4) else -1)
def chi8(n):
    r = n % 8
    return 0 if r%2==0 else (1 if r in (1,7) else -1)
N = 200000
s5 = mp.nsum(lambda n: chi5(int(n))/mp.mpf(n)**3, [1, N]) 
s8 = mp.nsum(lambda n: chi8(int(n))/mp.mpf(n)**3, [1, N])
print("partial sums:", mp.nstr(s5, 30), mp.nstr(s8, 30))
print("diffs (should be ~1/2N^2 = 1.25e-11):", L5-s5, L8-s8)
# eq:chi3 control
b = mp.sqrt(2)-1; L = mp.log(1+mp.sqrt(2))
chi3b = (Li3(b)-Li3(-b))/2
lhs = 48*chi3b; rhs = 32*mp.sqrt(2)*L8 - 3*mp.pi**2*L + 4*L**3
print("eq:chi3 check lhs-rhs:", lhs-rhs)
# golden beta
phi = (1+mp.sqrt(5))/2
D = L3(1/phi) - L3(-1/phi)
print("D/L5 =", D/L5)
import sympy as sp
print("nsimplify:", sp.nsimplify(str(D/L5)))

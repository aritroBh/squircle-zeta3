# step_26.py : FINAL summary evidence table: constants to 45 digits + residual tables
import math, pickle
from mpmath import mp, mpf, sqrt, log, nstr, pi, zeta
mp.dps = 70
L2 = log(2); z2 = pi**2/6; z3 = zeta(3)

beta = L2/2 - 1
gamma_ = -(1+sqrt(2)) + L2/2 - z2/2 - L2**2/8
b3 = mpf('0.5') - L2/2
c3 = L2**2/4 - L2/2 + 1 + z2
d3 = -L2**3/24 + L2**2/8 - L2/2 - L2*z2/2 + z2/2 - z3/2 + 1 + 3*sqrt(2)/2
A_D = L2
B_D = L2**2/2 + z2/2 - L2/2
C_D = L2**3/3 + L2*z2 + z3/2 - L2**2/2 - z2/2 + L2
print("beta  = (log2)/2 - 1                    =", nstr(beta, 45))
print("gamma = -(1+√2)+(log2)/2-π²/12-log²2/8  =", nstr(gamma_, 45))
print("b3    = (1-log2)/2                      =", nstr(b3, 45))
print("c3    = log²2/4-(log2)/2+1+π²/6         =", nstr(c3, 45))
print("d3    = -log³2/24+log²2/8-(log2)/2-(π²/12)log2+π²/12-ζ(3)/2+1+3√2/2")
print("      =", nstr(d3, 45))
print("A_D = log2 =", nstr(A_D, 45))
print("B_D = log²2/2+π²/12-(log2)/2 =", nstr(B_D, 45))
print("C_D = log³2/3+(π²/6)log2+ζ(3)/2-log²2/2-π²/12+log2 =", nstr(C_D, 45))
print()
print("m̃_1 = 1-log2 =", nstr(1-L2, 45))
print("m̃_2 = 1-log2+log²2/2+π²/12 =", nstr(1-L2+L2**2/2+z2/2, 45))
print()
d = pickle.load(open("E_values_new.pkl","rb"))
def R(n):
    L = log(mpf(n))
    return -L**2/2 + beta*L + gamma_
def R3(n):
    L = log(mpf(n))
    return L**3/3 + b3*L**2 + c3*L + d3
print("%-9s %-60s %-26s %-26s %-12s" % ("n", "E(n) (50 digits)", "resid = n²(E-√2/n)−R(n)", "R3num−R3(n)", "resid/(log³n/n)"))
for k,v in d.items():
    n = mpf(k); E = mpf(v[0]); resid = mpf(v[1]); L = log(n)
    print("%-9s %-60s %-26s %-26s %-12s" % ("1e%d"%int(round(log(n)/log(10))), nstr(E, 52), nstr(resid, 24), nstr(n*resid-R3(n), 22), nstr(resid/(L**3/n), 10)))

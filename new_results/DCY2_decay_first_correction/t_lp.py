
from flint import ctx, arb
ctx.dps = 40
from dcy2_cell import *

def dN_of(N):
    N = arb(N); kap = kap_of(N)
    return 2/(N*kap**2) + 1/(4*N)

def LamPrime(N):
    """d Lambda(-Delta(N))/dN, algebraically cancelled form."""
    N = arb(N)
    y = y_of(N); D = Delta_of(N); kap = kap_of(N)
    t = -D
    u = 1+y*t
    Q = u**2 + SQ2*y*u + y**2
    Q0 = 1+SQ2*y+y**2
    dN = dN_of(N)
    yN = -y/(2*N)
    tN = t*dN
    uN = yN*t + y*tN
    A0 = (u*u*Q0/Q).log()/2
    A0p = (SQ2*y*(u+SQ2*y)*uN/u - (SQ2*u+2*y)*yN)/(2*Q)
    D1p = -uN/u + bR*D*dN
    return A0 + N*A0p + D1p

for N in [arb(101000, 1000), arb(10100, 100), arb(2020, 20), arb(635, 7), arb(404, 4)]:
    v = LamPrime(N)
    print(f"cell mid={v.mid()}, rad={v.rad()}")
# compare with truth: FD
for N0 in [101000, 10100, 2020, 635, 404]:
    h = N0/1000 if N0 > 1000 else 1
    fd = (LogM_of(N0+h,-Delta_of(N0+h)) - LogM_of(N0-h,-Delta_of(N0-h)))/(2*h)
    print(f"N={N0}: FD={fd.str(10)}  LamPrime={LamPrime(arb(N0)).str(10)}")

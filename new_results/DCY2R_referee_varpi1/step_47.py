import sys
sys.path.insert(0, '../DCY2_decay_first_correction')
from flint import arb, acb, ctx
ctx.dps = 45
from step_07_partb import y_of, ell_of, H_of
PI = arb.pi(); SQ2 = arb(2).sqrt()
OMEGA = acb(1,1)/SQ2

def Hders(N):
    """H, H', H'', H''' at s0 via polygamma (same machinery as their certified momentsum)."""
    y = y_of(N); s0 = OMEGA/y
    H0 = H_of(N)
    Ls = []
    for j in range(1, 6):
        Ls.append( arb(2)**(-(j)) * ((s0/2).digamma() if j==1 else (s0/2).polygamma(j-1))
                   - arb(2)**(-(j)) * ((s0/2+arb(3)/2).digamma() if j==1 else (s0/2+arb(3)/2).polygamma(j-1)) )
    Hd = [H0]
    for m in range(1, 4):
        Hd.append(sum(arb(arb(m-1).fac()/(arb(k).fac()*arb(m-1-k).fac()))*Ls[k]*Hd[m-1-k] for k in range(m)))
    return s0, Hd

for N in [400, 3000, 100000]:
    y = y_of(N)
    s0, Hd = Hders(N)
    H0 = Hd[0]
    # g(s0) = H0 / (sqrt(pi/2) s0^{-3/2})
    g = H0/((PI/2).sqrt()*s0**(-arb(3)/2))
    r_g = g - 1 + 3/(4*s0)
    CG = abs(r_g)/y**2
    h1 = OMEGA*Hd[1]/H0
    h2 = OMEGA**2*Hd[2]/(2*H0)
    h3 = OMEGA**3*Hd[3]/(6*H0)
    r1 = abs(h1)/ (arb("1.5")*y)           # <= 1+2y ?
    r2 = abs(h2)/ (arb("1.875")*y**2)
    r3 = abs(h3)/ (arb("2.1875")*y**3)
    print(f"N={N:7d}: C_G = {CG.str(10)} (claim <=1)")
    print(f"          |h1|/1.5y = {r1.str(8)}, |h2|/1.875y^2 = {r2.str(8)}, |h3|/2.1875y^3 = {r3.str(8)} (claim <= 1+2y = {float(1+2*y):.5f})")

# c2 >= 0.9 bR y certified: c2/(bR y) as a function of y on interval [0, 0.063]
bR = PI/(2*SQ2)
def c2_ratio(y):
    # N = pi/(2 y^2); c2 = (N-1)y^2/2 - N y^2 (1+sqrt2 y)/(2 (1+sqrt2 y+y^2)^2)
    N = PI/(2*y**2)
    return ((N-1)*y**2/2 - N*y**2*(1+SQ2*y)/(2*(1+SQ2*y+y**2)**2))/(bR*y)
yv = arb(arb("0.0627")/2, arb("0.0627")/2 + arb("1e-4"))  # interval covering [0,0.0628]
# y appears as 1/y^2 -> singular at 0; instead check on [0.001, 0.0628] and note c2/(bR y) -> 1 as y->0
lo = c2_ratio(arb("0.0005")); hi = c2_ratio(arb("0.0628"))
print(f"c2/(bR y) at y=0.0005: {lo.str(10)}, at y=0.0628: {hi.str(10)}")
# interval evaluation over the whole range (y=0.0005..0.0628), N large => use monotonicity: check derivative sign numerically on grid
import math
worst = arb(1)
for i in range(200):
    yy = arb("0.0005") + (arb("0.0628")-arb("0.0005"))*i/199
    v = c2_ratio(yy)
    if v.lower() < worst: worst = v.lower()
print("min c2/(bR y) on grid:", worst.str(10), " (claim >= 0.9)")

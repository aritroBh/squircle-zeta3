# Referee check of thm:xnyn: area of z = x^n y^n over [0,1]^2. Own code, direct from the definition.
# A(n)-1 = iint (sqrt(1+G^2)-1) dx dy, G^2 = n^2 (xy)^{2n-2} (x^2+y^2).
# Coordinates x=e^{-s/n}, y=e^{-t/n}; sigma=s+t, delta=s-t (dx dy = n^-2 e^{-sigma/n} ds dt, ds dt = dsigma ddelta/2),
# symmetric in delta. No decomposition, no Mellin transforms.
from mpmath import mp, mpf, quad, sqrt, log, exp, pi, zeta, nstr, cosh, inf, fabs
import sys
mp.dps = 45
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s); sys.stdout.flush()

def E(n):
    n = mpf(n); e = 1 / n
    def inner(sig):
        # integrand in delta on [0, sig] (times 2 for symmetry, times 1/2 Jacobian -> 1)
        base = n * exp(-sig * (1 - e))
        def f(d):
            # e^{-2s/n}+e^{-2t/n} = 2 e^{-sig/n} cosh(d/n)
            G = base * sqrt(2 * exp(-sig * e) * cosh(d * e))
            G2 = G * G
            return G2 / (1 + sqrt(1 + G2))
        return quad(f, [0, sig]) * exp(-sig * e)
    s0 = log(n * sqrt(2))
    pts = [mpf(0)] + [s0 + k for k in range(-6, 7) if s0 + k > 0] + [s0 + 12, s0 + 25, s0 + 50, s0 + 90, s0 + 160]
    pts = sorted(set(pts))
    return quad(inner, pts) / n ** 2

lam = log(2); z3 = zeta(3)
b1 = lam / 2 - 1
b0 = -(1 + sqrt(2)) + lam / 2 - pi ** 2 / 12 - lam ** 2 / 8
b3 = (1 - lam) / 2
c3 = 1 - lam / 2 + lam ** 2 / 4 + pi ** 2 / 6
d3 = 1 + 3 * sqrt(2) / 2 - lam / 2 + lam ** 2 / 8 - lam ** 3 / 24 + pi ** 2 / 12 - pi ** 2 * lam / 12 - z3 / 2
pr("printed numerics: beta1 %s beta0 %s d3 %s" % (nstr(b1, 12), nstr(b0, 12), nstr(d3, 12)))

def pred(n):
    n = mpf(n); l = log(n)
    R = -l ** 2 / 2 + b1 * l + b0
    R3 = l ** 3 / 3 + b3 * l ** 2 + c3 * l + d3
    return sqrt(2) / n + R / n ** 2 + R3 / n ** 3, R, R3

# precision self-check
mp.dps = 60; Ea = E(1000); mp.dps = 45; Eb = E(1000)
pr("self-check n=1000: |E(dps45)-E(dps60)| / E =", nstr(fabs(Ea - Eb) / Ea, 3))

pr("n, E(n), r4 = n^4 (E - pred)/(1+log n)^4, and n^4(E-pred)/log^4 n")
for n in ['2', '3', '5', '10', '30', '100', '1e3', '1e4', '1e5', '1e6', '1e8', '1e10']:
    En = E(n); p, R, R3 = pred(n); nn = mpf(n)
    r = (En - p) * nn ** 4
    pr("  n=%s E=%s  r/(1+log n)^4=%s  r/log^4n=%s" % (n, nstr(En, 30), nstr(r / (1 + log(nn)) ** 4, 8), nstr(r / log(nn) ** 4, 8)))
    if n == '1e5':
        pr("   D(1e5)=n^2(E-sqrt2/n)=%s ; R+R3/n=%s" % (nstr(nn ** 2 * (En - sqrt(2) / nn), 10), nstr(R + R3 / nn, 10)))

# sensitivity: d3 without -zeta(3)/2
pr("dropping -zeta(3)/2 moves d3 by", nstr(z3 / 2, 6))
open(__file__.replace('.py', '.out'), 'w').write("\n".join(out) + "\n")

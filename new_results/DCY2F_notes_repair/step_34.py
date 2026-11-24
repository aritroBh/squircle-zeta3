import sys, time
sys.path.insert(0, "<project>/research/08_theory/DCY2F_notes_repair")
from flint import arb, acb, ctx
ctx.dps = 35
PI = arb.pi(); SQ2 = arb(2).sqrt(); W = acb(1,1)/SQ2
a_rate = PI/(2*SQ2)
def main_arb_fast(N):
    p = PI
    def f(r, analytic):
        s = r*W
        return ((acb(0,1)*p*s/2).exp() * (s/(1+s))**arb(N-2) / (1+s)**2
                * p.sqrt()/2 * (s/2+1).gamma() * (s/2+arb(3)/2).rgamma() * W)
    value = acb(0)
    pts = [0, 1, 4, 16, 64, 200]
    for lo, hi in zip(pts, pts[1:]):
        value += acb.integral(f, lo, hi, abs_tol=arb("1e-36"), rel_tol=arb("1e-32"), eval_limit=40000)
    R = arb(200)
    tail = 2*SQ2*(-a_rate*R).exp()/(p*a_rate*R**2)
    return 2/p*value.imag + arb(0, tail.upper())
for N in [100, 158, 250, 400]:
    t0=time.monotonic()
    M = main_arb_fast(N)
    dt = time.monotonic()-t0
    print(f"N={N}: M={M.str(12)}   ({dt:.2f}s)")

import sys, time
sys.path.insert(0, '../DCY2_decay_first_correction')
from flint import arb, acb, ctx
PI = arb.pi(); SQ2 = arb(2).sqrt()
W = acb(1,1)/SQ2
varpi1 = acb( (24+12*PI+PI**2)/(48*PI.sqrt()), (24-12*PI+PI**2)/(48*PI.sqrt()) )
a_rate = PI/(2*SQ2)

def main_fast(N):
    p = PI
    def f(r, analytic):
        s = r*W
        return ( (acb(0,1)*p*s/2).exp() * (s/(1+s))**arb(N-2) / (1+s)**2
                 * p.sqrt()/2 * (s/2+1).gamma() * (s/2+arb(3)/2).rgamma() * W )
    value = acb(0)
    pts = [0, 1, 4, 16, 64, 200]
    for lo, hi in zip(pts, pts[1:]):
        piece = acb.integral(f, lo, hi, abs_tol=arb("1e-36"), rel_tol=arb("1e-32"), eval_limit=30000)
        assert piece.is_finite()
        value += piece
    R = arb(200)
    tail = 2*SQ2*(-a_rate*R).exp()/(p*a_rate*R**2)
    return 2/p*value.imag + arb(0, tail.upper())

ctx.dps = 35
t0 = time.monotonic()
worst = arb(0)
for k in range(99, 400):
    N = arb(k+1)
    th = (PI*N).sqrt()
    M = main_fast(N)
    X = (PI.sqrt()/N)*(-th).exp()*((th).cos()*(1+varpi1.real/N.sqrt()) - (th).sin()*(varpi1.imag/N.sqrt()))
    res = abs(M + X) * N**2 * th.exp()
    if res.upper() > worst: worst = res.upper()
    assert res < arb(10)**6*N**(-2)*(-th).exp()
print(f"dps=35: 301 quadratures, max |M_main+X|N^2 e^th = {worst.str(8)}, elapsed {time.monotonic()-t0:.1f}s")

import sys, time
sys.path.insert(0, '../DCY2_decay_first_correction')
from flint import arb, acb, ctx
PI = arb.pi(); SQ2 = arb(2).sqrt()
W = acb(1,1)/SQ2
varpi1 = acb( (24+12*PI+PI**2)/(48*PI.sqrt()), (24-12*PI+PI**2)/(48*PI.sqrt()) )
a_rate = PI/(2*SQ2)
ctx.dps = 35

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
        value += piece
    R = arb(200)
    tail = 2*SQ2*(-a_rate*R).exp()/(p*a_rate*R**2)
    return 2/p*value.imag + arb(0, tail.upper())

t0=time.monotonic()
bad = []
worst = arb(0); wN=None
for k in range(99, 400):
    N = arb(k+1)
    M = main_fast(N)
    th = (PI*N).sqrt()
    X = (PI.sqrt()/N)*(-th).exp()*((th).cos()*(1+varpi1.real/N.sqrt()) - (th).sin()*(varpi1.imag/N.sqrt()))
    res = abs(M + X) * N**2 * th.exp()
    if res.upper() > worst: worst = res.upper(); wN = k+1
    if not (res < arb(10)**6*N**(-2)*(-th).exp()):
        bad.append((k+1, res.str(8)))
print("failures:", bad[:10])
print(f"max res = {worst.str(8)} at N={wN}, elapsed {time.monotonic()-t0:.1f}s")

"""step_06_finite.py — certified M_main at each integer N in [100,400] (arb quadrature).

M_main(N) = (2/pi) Im int_0^inf H(omega x) e^{i pi omega x /2} (omega x/(1+omega x))^N omega dx / (omega x),
H(s) = (1/2) Beta(s/2, 3/2) = sqrt(pi)/4 Gamma(s/2)/Gamma(3/2+s/2).
Segments [0,1,4,16,64,200] + proved tail  2 sqrt2 e^{-aR}/(pi a R^2), a = pi/(2 sqrt2), R = 200
(CX6 proof.md eq. (20)).
Then: residual  |M_main + (sqrt pi/N) e^{-theta} Re[e^{i theta}(1+varpi1/sqrt N)]|
       <= C2 N^{-2} e^{-theta}   checked with arb (outward).
Also cross-check vs mk_deg9.tsv (certified M_k balls).
"""
from flint import arb, acb, ctx
from pathlib import Path

PI = arb.pi(); SQ2 = arb(2).sqrt()
W = acb(1,1)/SQ2
a_rate = PI/(2*SQ2)
varpi1 = acb( (24+12*PI+PI**2)/(48*PI.sqrt()), (24-12*PI+PI**2)/(48*PI.sqrt()) )

def main_arb(N):
    p = PI
    def f(r, analytic):
        s = r*W
        return ( (acb(0,1)*p*s/2).exp() * (s/(1+s))**arb(N-2) / (1+s)**2
                 * p.sqrt()/2 * (s/2+1).gamma() * (s/2+arb(3)/2).rgamma() * W )
    value = acb(0)
    pts = [0, 1, 4, 16, 64, 200]
    for lo, hi in zip(pts, pts[1:]):
        piece = acb.integral(f, lo, hi, abs_tol=arb("1e-46"), rel_tol=arb("1e-40"), eval_limit=80000)
        assert piece.is_finite()
        value += piece
    R = arb(200)
    tail = 2*SQ2*(-a_rate*R).exp()/(p*a_rate*R**2)
    return 2/p*value.imag + arb(0, tail.upper())

def residual_check(N):
    N = arb(N)
    th = (PI*N).sqrt()
    M = main_arb(N)
    T1 = -(PI.sqrt()/N)*(-th).exp()*( (th).cos()*(1+varpi1.real/N.sqrt()) - (th).sin()*(varpi1.imag/N.sqrt()) )
    res = abs(M - T1)
    bound = arb(10)**6 * N**(-2) * (-th).exp()
    return M, res, bound

if __name__ == "__main__":
    ctx.dps = 45
    ROOT = Path("<project>")
    rows = {}
    for line in (ROOT/"research/08_theory/D13f_sweep/mk_deg9.tsv").read_text().splitlines():
        f = line.split("	")
        if f[0].isdigit():
            rows[int(f[0])] = arb(f[1])
    import time
    t0 = time.monotonic()
    worst = arb(0); worstN = None
    for k in range(99, 400):   # N = k+1 in [100, 400]
        N = arb(k+1)
        M, res, bound = residual_check(N)
        # table cross-check for the full M_k statement (sharp): |M_k - T1| <= res + |M_k - M_main|
        Mk = rows[k]
        dev = abs(Mk - M)
        res_full = abs(Mk - (-(PI.sqrt()/N)*(-(PI*N).sqrt()).exp()*(
            (PI*N).sqrt().cos()*(1+varpi1.real/N.sqrt()) - (PI*N).sqrt().sin()*(varpi1.imag/N.sqrt()))))
        r = res/bound
        if r.upper() > worst:
            worst = r.upper(); worstN = N
        if k % 50 == 0:
            print(f"N={N}: M_main={M.str(12)}")
            print(f"   res={res.str(6)}  bound={bound.str(6)}  res/bound={r.str(6)}  |Mk-Mmain|={dev.str(6)}")
        assert res < bound, (N, res, bound)
    print(f"ALL PASS: max res/bound = {worst.str(10)} at N={worstN}; elapsed {time.monotonic()-t0:.1f}s")

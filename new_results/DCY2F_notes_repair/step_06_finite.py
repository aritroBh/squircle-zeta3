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

def X_of(N):
    # X = (sqrt pi/N) e^{-theta} Re[e^{i theta}(1+varpi1/sqrt N)]  -- the theorem's main-term term.
    # (In the old notation of verify_dcy2.py/step_06_finite.py, T1 = -X; note the SIGN:
    #  the theorem needs |M_main + X|, NOT |M_main - X| ~ 2|X|.)
    N = arb(N); th = (PI*N).sqrt()
    return (PI.sqrt()/N)*(-th).exp()*( (th).cos()*(1+varpi1.real/N.sqrt()) - (th).sin()*(varpi1.imag/N.sqrt()) )

def residual_check(N):
    N = arb(N)
    th = (PI*N).sqrt()
    M = main_arb(N)
    X = X_of(N)
    res = abs(M + X)   # R1 (referee DCY2R): theorem quantity is |M_main + X|
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
    worst = arb(0); worstN = None; worstF = arb(0); worstNF = None
    for k in range(99, 400):   # N = k+1 in [100, 400]
        N = arb(k+1)
        M, res, bound = residual_check(N)          # res = |M_main + X|  (R1)
        # full-M_k sharp residual from the certified table: |M_k + X|
        dev = abs(rows[k] - M)
        res_full = abs(rows[k] + X_of(N))
        r = res/bound
        if r.upper() > worst:
            worst = r.upper(); worstN = N
        if res_full > worstF:
            worstF = res_full.upper(); worstNF = N
        if k % 50 == 0:
            print(f"N={N}: M_main={M.str(12)}")
            print(f"   |M_main+X|={res.str(6)}  bound={bound.str(6)}  res/bound={r.str(6)}  |Mk-Mmain|={dev.str(6)}")
        assert res < bound, (N, res, bound)
    scale = lambda NN: arb(NN)**2*((PI*arb(NN)).sqrt()).exp()
    print(f"ALL PASS: max |M_main+X|/(C2 N^-2 e^-th) = {worst.str(10)} at N={worstN}")
    print(f"   i.e. max |M_main+X| N^2 e^th = {(worst*arb(10)**6).str(10)}  (constant needed, ~1.75)")
    print(f"   full M_k: max |M_k+X| N^2 e^th = {worstF.str(10)} at N={worstNF}")
    print(f"   elapsed {time.monotonic()-t0:.1f}s")

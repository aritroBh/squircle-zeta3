import sys, time
sys.path.insert(0, '../DCY2_decay_first_correction')
from flint import arb, acb, ctx
ctx.dps = 50
PI = arb.pi(); SQ2 = arb(2).sqrt()
W_ = acb(1,1)/SQ2
from step_07_partb import momentsum_of, H_of, y_of

def J_arb(N, segs):
    p = PI
    def f(r, analytic):
        s = r*W_
        return ( (acb(0,1)*p*s/2).exp() * (s/(1+s))**arb(N-2) / (1+s)**2
                 * p.sqrt()/2 * (s/2+1).gamma() * (s/2+arb(3)/2).rgamma() * W_ )
    value = acb(0)
    for lo, hi in zip(segs, segs[1:]):
        piece = acb.integral(f, lo, hi, abs_tol=arb(10)**(-45), rel_tol=arb("1e-40"), eval_limit=120000)
        assert piece.is_finite()
        value += piece
    R = arb(segs[-1])
    a_rate = p/(2*SQ2)
    tail = 2*SQ2*(-a_rate*R).exp()/(p*a_rate*R**2)
    return value + acb(arb(0, tail.upper()))

def S_of_exp(N):
    y = y_of(N)
    # e^{phi(s0)} = -i y e^{-th+i th} e^{S(y)}; we need model = (2/pi) e^{phi} H0 omega MS
    # = -i y e^{-th+i th} (2/pi) e^{S} H0 omega MS ; compute via exact phi instead:
    s0 = W_/y
    N_ = arb(N)
    phi = (N_-1)*acb(s0).log() - N_*(1+s0).log() + acb(0,1)*PI*s0/2
    return phi

t0=time.monotonic()
for N, segs in [(10000, [0,1,4,16,64,200]), (100000, [0,1,4,16,64,200,500,1000])]:
    J = (2/PI)*J_arb(N, segs)
    phi = S_of_exp(N)
    model = (2/PI)*phi.exp()*H_of(N)*W_*momentsum_of(N)
    N_ = arb(N); th = (PI*N_).sqrt()
    d = abs(J - model)*N_**2*th.exp()
    print(f"N={N}: |J-model|N^2 e^th (CERTIFIED arb) = {d.str(10)}   [J={J.str(6)}]")
print("elapsed", time.monotonic()-t0)

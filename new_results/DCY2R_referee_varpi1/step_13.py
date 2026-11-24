import sys
sys.path.insert(0, '../DCY2_decay_first_correction')
from flint import arb, acb, ctx
ctx.dps = 45
PI = arb.pi(); SQ2 = arb(2).sqrt()
W = acb(1,1)/SQ2
varpi1 = acb( (24+12*PI+PI**2)/(48*PI.sqrt()), (24-12*PI+PI**2)/(48*PI.sqrt()) )

def Itilde(N):
    p = PI
    def f(r, analytic):
        s = r*W
        return ( (acb(0,1)*p*s/2).exp() * (s/(1+s))**arb(N-2) / (1+s)**2
                 * p.sqrt()/2 * (s/2+1).gamma() * (s/2+arb(3)/2).rgamma() * W )
    value = acb(0)
    pts = [0, 1, 4, 16, 64, 200]
    for lo, hi in zip(pts, pts[1:]):
        value += acb.integral(f, lo, hi, abs_tol=arb("1e-46"), rel_tol=arb("1e-40"), eval_limit=60000)
    return value   # tail negligible for these N? check separately

for N in (200, 400):
    N_ = arb(N); th = (PI*N_).sqrt()
    I = Itilde(N)
    T1 = -acb(0,1)*(PI.sqrt()/N_)*((-th+i*th) if False else 0)
    T1 = acb(0,-1)*(PI.sqrt()/N_)* (acb(th,th)).exp() * (1+varpi1/N_.sqrt())
    # T1 = -i (sqrt pi /N) e^{-th + i th} (1+varpi1/sqrt N):  e^{-th+ith} = exp(acb(-th,th))
    T1 = acb(0,-1)*(PI.sqrt()/N_)* acb(-th,th).exp() * (1+varpi1/N_.sqrt())
    print(f"N={N}")
    print("  Itilde        =", I.str(12))
    print("  T1            =", T1.str(12))
    print("  (pi/2) T1     =", (T1*PI/2).str(12))
    print("  |I-T1|/(|T1|) =", (abs(I-T1)/abs(T1)).str(6))
    print("  |I-(pi/2)T1|/(|T1|) =", (abs(I-T1*PI/2)/abs(T1)).str(6))

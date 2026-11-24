import sys, time
sys.path.insert(0, '../DCY2_decay_first_correction')
from flint import arb, acb, ctx
ctx.dps = 45
t0 = time.monotonic()
from step_06_finite import main_arb, varpi1
PI = arb.pi()
for N in (100, 158, 200, 381, 400):
    N_ = arb(N)
    th = (PI*N_).sqrt()
    M = main_arb(N_)
    X = (PI.sqrt()/N_)*(-th).exp()*((th).cos()*(1+varpi1.real/N_.sqrt()) - (th).sin()*(varpi1.imag/N_.sqrt()))
    # theorem quantity: |M_main + X| ; DCY2 measured: |M_main - X|
    res_thm  = abs(M + X) * N_**2 * th.exp()
    res_dcy2 = abs(M - X) * N_**2 * th.exp()
    print(f"N={N}: M_main={M.str(8)}")
    print(f"   theorem  |M+X| N^2 e^th = {res_thm.str(8)}")
    print(f"   DCY2's   |M-X| N^2 e^th = {res_dcy2.str(8)}   (2*sqrt(pi)*N = {float(2*PI.sqrt()*N_):.2f})")
print(f"elapsed {time.monotonic()-t0:.1f}s")

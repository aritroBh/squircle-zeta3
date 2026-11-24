import sys, time
sys.path.insert(0, '../DCY2_decay_first_correction')
from flint import arb, acb, ctx
from pathlib import Path
ctx.dps = 45
t0 = time.monotonic()
from step_06_finite import main_arb, varpi1
PI = arb.pi()
C2 = arb(10)**6

ROOT = Path("<project>")
rows = {}
for line in (ROOT/"research/08_theory/D13f_sweep/mk_deg9.tsv").read_text().splitlines():
    f = line.split("\t")
    if f[0].isdigit():
        rows[int(f[0])] = arb(f[1])

worst_main = arb(0); wN1=None
worst_full = arb(0); wN2=None
worst_ratio = arb(0)
for k in range(99, 400):
    N = arb(k+1)
    th = (PI*N).sqrt()
    X = (PI.sqrt()/N)*(-th).exp()*((th).cos()*(1+varpi1.real/N.sqrt()) - (th).sin()*(varpi1.imag/N.sqrt()))
    M = main_arb(N)
    res_main = abs(M + X) * N**2 * th.exp()          # theorem quantity for M_main
    Mk = rows[k]
    res_full = abs(Mk + X) * N**2 * th.exp()          # full M_k from certified table
    if res_main.upper() > worst_main: worst_main = res_main.upper(); wN1 = k+1
    if res_full.upper() > worst_full: worst_full = res_full.upper(); wN2 = k+1
    r = res_full/C2
    if r.upper() > worst_ratio: worst_ratio = r.upper()
print(f"[REPAIR, correct sign] N in [100,400]:")
print(f"  max |M_main + X| N^2 e^theta = {worst_main.str(10)} at N={wN1}")
print(f"  max |M_k    + X| N^2 e^theta = {worst_full.str(10)} at N={wN2}")
print(f"  max res_full/bound (C2=1e6)  = {worst_ratio.str(10)}")
print(f"  elapsed {time.monotonic()-t0:.1f}s")
open('progress.log','a').write(f"step: certified finite-range repair, max full-Mk residual {worst_full.str(6)} at N={wN2}\n")

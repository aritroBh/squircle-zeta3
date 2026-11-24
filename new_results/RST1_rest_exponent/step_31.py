import os
os.chdir('<project>/research/08_theory/RST1_rest_exponent')
step = r'''"""RST1 step_16: certified M_main quadrature timing test + sanity for a few k."""
from flint import arb, acb, ctx
from pathlib import Path
import time
ctx.dps = 45
PI = arb.pi(); SQ2 = arb(2).sqrt()
W = acb(1,1)/SQ2
a_rate = PI/(2*SQ2)
sq2pi = (2*PI).sqrt()

def main_arb(k):
    N = k+1
    p = PI
    def f(r, analytic):
        s = r*W
        return ( (acb(0,1)*p*s/2).exp() * (s/(1+s))**arb(N-2) / (1+s)**2
                 * p.sqrt()/2 * (s/2+1).gamma() * (s/2+arb(3)/2).rgamma() * W )
    value = acb(0)
    pts = [0, 1, 4, 16, 64, 200]
    for lo, hi in zip(pts, pts[1:]):
        piece = acb.integral(f, lo, hi, abs_tol=arb("1e-60"), rel_tol=arb("1e-55"), eval_limit=80000)
        assert piece.is_finite()
        value += piece
    R = arb(200)
    tail = 2*SQ2*(-a_rate*R).exp()/(p*a_rate*R**2)
    return 2/p*value.imag + arb(0, tail.upper())

ROOT = Path("<project>")
rows = {}
for line in (ROOT/"research/08_theory/D13f_sweep/mk_deg9.tsv").read_text().splitlines():
    f_ = line.split("\t")
    if f_[0].isdigit(): rows[int(f_[0])] = arb(f_[1])

C = arb(1030)
t0=time.monotonic()
for k in [99, 400]:
    N = arb(k+1)
    Mm = main_arb(k)
    rest = abs(rows[k]-Mm)
    tgt = C*N**arb('-0.5')*(-sq2pi*N.sqrt()).exp()
    th2 = sq2pi*N.sqrt()
    lead = PI.sqrt()/N*(-th2).exp()
    print(f"k={k}: M_main={Mm.str(15)}", flush=True)
    print(f"   |rest|={rest.str(12)}  C*target={tgt.str(12)}  ratio rest/(C tgt)={(rest/tgt).str(6)}", flush=True)
    print(f"   |rest|/target={(rest/(tgt/C)).str(6)}  |rest|/lead={(rest/lead).str(6)}  |cos th2|={abs(th2.cos()).str(6)}", flush=True)
print(f"elapsed for 2 k: {time.monotonic()-t0:.2f}s")
'''
open('step_16.py','w').write(step)
import subprocess
r = subprocess.run(['bash','-c','cd <project>/research/08_theory/RST1_rest_exponent && nice <project>/.venv/bin/python step_16.py'], capture_output=True, text=True, timeout=2400)
print(r.stdout)
print(r.stderr[-1500:] if r.stderr else '')

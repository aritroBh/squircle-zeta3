"""RST1 step_17: FULL certified sanity table, k=99..400 (N in [100,401]).
Checks |M_k - M_main| <= 1030 N^{-1/2} e^{-sqrt(2 pi N)}  (arb, outward).
Prints also |rest|/target and |rest|/lead (predicted-ratio evidence).
Checkpoint: sanity_progress.json (resume-safe).
"""
from flint import arb, acb, ctx
from pathlib import Path
import time, json, os
ctx.dps = 45
PI = arb.pi(); SQ2 = arb(2).sqrt()
W = acb(1,1)/SQ2
a_rate = PI/(2*SQ2); sq2pi=(2*PI).sqrt()

def main_arb(k):
    N = k+1; p = PI
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
ck = Path("sanity_progress.json")
done = json.loads(ck.read_text()) if ck.exists() else {}
t0 = time.monotonic()
maxratio = arb(0); maxk = None
for k in range(99, 401):
    if str(k) in done:
        r_ = arb(done[str(k)]["rest_mid"]), arb(done[str(k)]["rest_rad"])
        rest = arb(done[str(k)]["rest_mid"]) + arb(0, done[str(k)]["rest_rad"])
    else:
        Mm = main_arb(k)
        rest = abs(rows[k]-Mm)
        done[str(k)] = {"rest_mid": rest.mid().str(30), "rest_rad": rest.rad().str(30)}
        if k % 10 == 0: ck.write_text(json.dumps(done))
    N = arb(k+1)
    tgt = C*N**arb('-0.5')*(-sq2pi*N.sqrt()).exp()
    assert rest < tgt, (k, rest, tgt)
    rr = (rest/(tgt/C)).upper()
    if rr > maxratio: maxratio = rr; maxk = k
    if k % 25 == 0:
        th2 = sq2pi*N.sqrt()
        print(f"k={k}: |rest|={rest.str(10)}  |rest|/target={(rest/(tgt/C)).str(6)}  |rest|/lead={(rest/(PI.sqrt()/N*(-th2).exp())).str(6)}  cos_th2={th2.cos().str(6)}", flush=True)
ck.write_text(json.dumps(done))
print(f"ALL PASS k=99..400: |M_k-M_main| < 1030 N^-1/2 e^-sqrt(2piN).", flush=True)
print(f"max |rest|/target = {maxratio.str(8)} at k={maxk};  elapsed {time.monotonic()-t0:.1f}s", flush=True)

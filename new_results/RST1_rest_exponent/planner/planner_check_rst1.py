# Planner check of RST1 (independent of RST1's own code): rest = exact M_k (D13f certified balls) minus M_main from DCY2's
# refereed quadrature (step_06_finite.main_arb). Checks the claimed bound 1030 N^{-1/2} e^{-sqrt(2 pi N)} and the claimed
# leading law rest ~ (sqrt(pi)/N) e^{-sqrt(2 pi N)} cos sqrt(2 pi N) (envelope ratio A_k should track cos, |A_k| = O(1)).
import sys, time
from pathlib import Path
ROOT = Path("<project>")
sys.path.insert(0, str(ROOT/"research/08_theory/DCY2_decay_first_correction"))
from flint import arb, ctx
import step_06_finite as d
ctx.dps = 45
rows = {}
for line in (ROOT/"research/08_theory/D13f_sweep/mk_deg9.tsv").read_text().splitlines():
    f = line.split("\t")
    if f[0].isdigit(): rows[int(f[0])] = arb(f[1])
PI = arb.pi(); ok = True; t0 = time.monotonic()
for k in [99, 106, 130, 150, 175, 200, 250, 300, 350, 399]:
    N = arb(k+1); M = d.main_arb(k+1); rest = rows[k] - M
    th2 = (2*PI*N).sqrt()
    bound = 1030 * N**arb(-0.5) * (-th2).exp()
    A = rest * N * th2.exp() / PI.sqrt()
    alt = rest * (2*N.sqrt()).exp()              # if the true exponent were 2 (CX6), this would stay O(1)
    good = bool(abs(rest) < bound); ok &= good
    print(f"k={k}: rest={rest.str(8)}  rest/bound={(abs(rest)/bound).str(4)}  A_k={A.str(5)}  cos(th2)={th2.cos().str(5)}  rest*e^(2sqrtN)={alt.str(3)}  {'PASS' if good else 'FAIL'}")
print("ALL PASS" if ok else "SOME FAIL", f"({time.monotonic()-t0:.0f}s)")

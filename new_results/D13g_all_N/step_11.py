
import os
from flint import arb, ctx
ctx.dps = 60
BASE='research/08_theory'
lines = open(os.path.join(BASE,'D13f_sweep','mk_deg9.tsv')).read().strip().split('\n')
def parbarb(s):  # "[x +/- y]" -> arb mid as arb
    mid = s.strip()[1:].split(' ')[0]
    return arb(mid)
pi = arb.pi()
truth = {}
for ln in lines:
    p = ln.split('\t')
    k = int(p[0]); N = k+1
    Mk = parbarb(p[1])
    th = (pi*N).sqrt()
    lead = -(pi.sqrt()/N)*(-th).exp()*(th.cos())   # + (sqpi/N) e^{-th} cos th term added
    # theorem: |Mk + (sqpi/N) e^{-th} cos th| <= C N^{-3/2} e^{-th}
    resid = Mk + (pi.sqrt()/N)*(-th).exp()*(th.cos())
    Ctrue = abs(resid)*arb(N)**arb("1.5")*(th).exp()
    truth[N] = arb(Ctrue.upper())
mx = max(truth[N] for N in range(100,401))
Nmx = max(range(100,401), key=lambda N: truth[N])
print("max true C over [100,400]:", mx, "at N =", Nmx)
for N in [100,150,200,300,400]:
    print(N, truth[N])

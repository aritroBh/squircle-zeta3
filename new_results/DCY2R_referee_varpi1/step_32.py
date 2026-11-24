import sys, time
sys.path.insert(0, '../DCY2_decay_first_correction')
from mpmath import mp, mpf, pi, sqrt, exp, log
mp.dps = 40
import dcy2_audit as A
import dcy2_budget as B
t0=time.monotonic()
B.P.kap0 = 2.48
# truth |It - Model| scaled by N^2 e^theta, vs certified LHS_of (arb) -> use mpmath mirror of LHS pieces via budget quant
for N in [100000, 300000, 1000000, 3000000, 10000000]:
    tq = A.true_quant(N)
    truth = abs(tq['It'] - tq['Model'])*N**2*exp(tq['theta'])
    qq = B.quant(N)
    protoLHS = qq['E_total']*N**2*exp(tq['theta'])
    print(f"N={N:9d}: truth|J-model|N^2e^th = {float(truth):12.4f}   prototype E_total scaled = {float(protoLHS):12.4f}   ratio = {float(truth/protoLHS):.4f}")
print("elapsed", time.monotonic()-t0)

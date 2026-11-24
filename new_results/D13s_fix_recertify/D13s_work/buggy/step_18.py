# step_15.py - Task 2: tighten [100,400] with ratio-1.05 cells via step_38
import math, time, json
from step_38 import certified_y38
ys = lambda N: math.sqrt(math.pi/(2*N))
t0 = time.time()
r = 1.05
cells = []
N0 = 100.0
while N0 < 400*(1-1e-12):
    N1 = min(N0*r, 400.0)
    cells.append((N0, N1))
    N0 = N1
print(f"{len(cells)} cells, ratios <= {r}")
out = []
supC = 0; worst = None
for (n1, n2) in cells:
    res = certified_y38(ys(n2), ys(n1), J=24, nsub=4)
    C = float(res['C'].upper())
    pieces = {k: float(v.upper()) for k,v in res.items() if k!='C'}
    out.append([n1, n2, C, pieces])
    if C > supC: supC = C; worst = (n1, n2)
    print(f"[{n1:7.2f},{n2:7.2f}] C={C:9.3f}  W1={pieces['W1']:7.2f} W2={pieces['W2']:6.2f} W3={pieces['W3']:7.2f} TM={pieces['TM']:6.2f} TR={pieces['TR']:5.2f} TLm={pieces['TLm']:6.2f} TLf={pieces['TLf']:6.2f} D={pieces['D']:5.2f}")
print(f"\nsup C on [100,400] = {supC:.4f} at {worst}   ({time.time()-t0:.1f}s)")
json.dump(out, open("tight_100_400.json","w"))
with open('progress.log','a') as f:
    f.write(f"turn4: [100,400] ratio-1.05 sweep: {len(cells)} cells, supC={supC:.4f} at {worst}.\n")

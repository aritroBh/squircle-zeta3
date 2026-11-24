# Deterministic adaptive driver around D13g's certified_y38 (no new mathematics): cover [400, 1e10] by N-cells,
# halving the log-ratio of a cell when an assert fails or C > CMAX. Prints every accepted cell and the sup.
import math, json, sys, time
RATIO = float(sys.argv[1]) if len(sys.argv)>1 else 1.05
sys.path.insert(0, '.')
from step_38 import certified_y38
from flint import arb
ys = lambda N: math.sqrt(math.pi/(2*N))
CMAX = float("inf"); out = []; sup = 0.0; fails = 0
stack = [(400.0, 1e10)]
cells = []
# initial ratio-1.5 grid
n = 400.0
while n < 1e10:
    m = min(n*RATIO, 1e10); cells.append((n, m)); n = m
t0 = time.time()
while cells:
    n1, n2 = cells.pop(0)
    try:
        r = certified_y38(ys(n2), ys(n1), nsub=4)
        C = float(r['C'].upper())
        ok = math.isfinite(C) and C <= CMAX
    except AssertionError as e:
        ok = False; C = float('nan')
    if ok:
        out.append((n1, n2, C)); sup = max(sup, C)
    else:
        fails += 1
        if n2/n1 < 1.0005: print(f'GIVE UP cell [{n1:.6g},{n2:.6g}] C={C}'); out.append((n1, n2, None)); continue
        mid = math.sqrt(n1*n2); cells[:0] = [(n1, mid), (mid, n2)]
print(f'cells accepted {sum(1 for o in out if o[2] is not None)}, gave up {sum(1 for o in out if o[2] is None)}, splits {fails}, time {time.time()-t0:.1f}s')
print(f'sup C over accepted cells = {sup:.6g}')
covered = all(o[2] is not None for o in out) and abs(min(o[0] for o in out)-400)<1e-9 and abs(max(o[1] for o in out)-1e10)<1
print('contiguous cover of [400,1e10] with all cells certified:', covered)
for o in out[::max(1,len(out)//12)]: print('  ', f'[{o[0]:.6g},{o[1]:.6g}]', o[2])
json.dump(out, open('sweep_cells_r'+str(RATIO)+'.json','w'))

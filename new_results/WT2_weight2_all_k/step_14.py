import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
import sympy as sp
from ibp_reduce import reduce_k

v = sp.Symbol('v', positive=True)
S2 = sp.sqrt(2); a = 1 + S2; b = S2 - 1
Rw2 = 1/(2*v) + 1/(2*(v - a)) - 1/(v + a) - 1/(2*(v + b))

loglines = []
def log(msg):
    print(msg, flush=True)
    loglines.append(msg)

for k in [2, 3, 4, 5]:
    t0 = time.time()
    res = reduce_k(k, verbose=False)
    dt = time.time() - t0
    log(f"k = {k}: keys {sorted(res.keys())}  [{dt:.2f}s]")
    for key, val in res.items():
        if sum(key) != 1:
            continue
        d = sp.simplify(sp.radsimp(sp.together(val - Rw2)))
        log(f"   weight-2 key {key}:  kernel - Rw2 = {d}")
    deg0 = {kk: vv for kk, vv in res.items() if sum(kk) == 0}
    log(f"   degree-0 (weight-1) pieces: {deg0}")

with open('progress.log', 'a') as f:
    f.write("\n".join(loglines) + "\n")

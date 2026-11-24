
from flint import ctx, arb
ctx.dps = 40
import dcy2_cert as C
for N in [400, 2000, 10000, 10**6]:
    r, q = C.ratio(N)
    print(f"N={N:9d} LHS/C2={r.str(10)}  LHS={q['LHS'].str(8)}")
print("--- cells")
for lo, hi in [(400, 520), (520, 676), (1000, 1300), (2000, 2600), (10000, 13000), (10**5, 1.3e5), (10**6, 1.3e6)]:
    N = arb((lo+hi)/2, (hi-lo)/2)
    r, q = C.ratio(N)
    print(f"[{lo}, {hi}]  ratio = {r.str(10)}")


from flint import ctx
ctx.dps = 40
import dcy2_cell as C
C2 = 10**6
for N in [400, 2000, 3000, 10000, 10**5, 10**6, 10**9, 10**12]:
    v = C.LHS_of(N)
    print(f"N={N:9d}  LHS={v.str(10)}   LHS/C2={float(v/C2):.5f}")


from flint import ctx, arb
ctx.dps = 40
import dcy2_cert as C
N = arb(460, 60)
q = C.quantities(N)
for k, v in q.items():
    try:
        print(k, "=", v.str(10) if hasattr(v, "str") else v)
    except Exception as e:
        print(k, "ERR", e)

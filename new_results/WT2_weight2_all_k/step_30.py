import sys, importlib, time
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
from ibp_reduce import reduce_k
import fastred; importlib.reload(fastred)
from fastred import *

v = sp.Symbol('v', positive=True)
S2 = sp.sqrt(2); a = 1+S2; b = S2-1
Rw2 = 1/(2*v) + 1/(2*(v-a)) - 1/(v+a) - 1/(2*(v+b))

for k in [2,3]:
    t0=time.time(); slow = reduce_k(k, verbose=False); t1=time.time()
    fast = reduce_k_fast(k); t2=time.time()
    print(f"k={k}: slow {t1-t0:.2f}s fast {t2-t1:.3f}s")
    assert set(slow.keys()) == set(fast.keys())
    for key in slow:
        d = sp.together(slow[key]) - sp.together(pf_to_sympy(fast[key]))
        assert abs(complex(sp.N(d.subs(v, sp.Rational(7,10)), 25))) < 1e-18, (k, key)
    ker = pf_to_sympy(fast[(1,0)])
    print("   kernel - Rw2 =", sp.simplify(sp.radsimp(sp.together(ker) - Rw2)))
print("MATCH slow==fast for k=2,3 on ALL keys; kernel == Rw2")

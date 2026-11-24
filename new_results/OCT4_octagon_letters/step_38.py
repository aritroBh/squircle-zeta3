import os, time; os.chdir('research/08_theory/OCT4_octagon_letters')
from oct4_common import log
import sympy as sp
t = sp.symbols('t', positive=True)
a, b = sp.symbols('a b')
# Try the definite integral of log(1+a t)^2 log(1+b t) — symbolic
t0=time.time()
try:
    prim = sp.integrate(sp.log(1+a*t)**2*sp.log(1+b*t), t)
    print("indefinite integral computed in", round(time.time()-t0,1),"s")
    print(sp.simplify(prim))
except Exception as e:
    print("FAIL", e, time.time()-t0)
log("step: sympy symbolic integration attempt for cross integral")

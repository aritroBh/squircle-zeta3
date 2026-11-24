import os, time; os.chdir('research/08_theory/OCT4_octagon_letters')
from oct4_common import log as flog
import sympy as sp
t=sp.symbols('t'); a,b=sp.symbols('a b', nonzero=True)
t0=time.time()
# U = int log(t) log(1+b t)/(1+a t) dt
expr = sp.log(t)*sp.log(1+b*t)/(1+a*t)
try:
    U = sp.integrate(expr, t)
    print("U indefinite computed", round(time.time()-t0,1),"s")
    print(sp.simplify(U))
except Exception as e:
    print("FAIL",repr(e),round(time.time()-t0,1))
flog("step: sympy integrate U")

import sympy as sp, mpmath, json, pickle, time
from mpmath import mp, mpf, pi, nstr, log, zeta, polylog
mp.dps = 110
u, g = sp.symbols('u g', positive=True)
ca = sp.symbols('ca')
fks = pickle.load(open('ckpt_corner_integrand.pkl','rb'))['fks']
import math
t0=time.time()
res = {}
for alpha_name, alpha_val in [('pi/2', mp.pi/2), ('3pi/4', 3*mp.pi/4), ('2pi/3', 2*mp.pi/3)]:
    cval = mp.cos(alpha_val)
    vals = []
    for k, fk in enumerate(fks):
        fkl = sp.lambdify((g, ca), fk, 'mpmath')
        fun = lambda x: fkl(x, cval)
        v = mpmath.quad(fun, [mpf(0), mpf('0.5'), mpf(1)])
        vals.append(nstr(v, 100))
        print(f"alpha={alpha_name} k={k}: {nstr(v,50)}   ({time.time()-t0:.0f}s)")
    res[alpha_name] = vals
json.dump(res, open('ckpt_Mk_cornerintegral.json','w'))
with open('progress.log','a') as f:
    f.write("step_32: M_k(alpha) k=0..3 at pi/2, 3pi/4, 2pi/3 by direct quadrature of the u-expanded corner integrand (110 dps)\n")

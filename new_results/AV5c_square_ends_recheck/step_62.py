import sympy as sp, mpmath, json, pickle, time
from mpmath import mp, mpf, pi, nstr
mp.dps = 300
u, g = sp.symbols('u g', positive=True)
ca = sp.symbols('ca')
fks = pickle.load(open('ckpt_corner_integrand.pkl','rb'))['fks']
cval = mp.cos(3*mp.pi/4/2)  # cos(3pi/8)
t0 = time.time()
vals = []
for k, fk in enumerate(fks):
    fkl = sp.lambdify((g, ca), fk, 'mpmath')
    fun = lambda x: fkl(x, cval)
    v = mpmath.quad(fun, [mpf(0), mpf('0.25'), mpf('0.5'), mpf('0.75'), mpf(1)])
    vals.append(nstr(v, 280))
    print(f"k={k}: {nstr(v,60)}   ({time.time()-t0:.0f}s)")
json.dump({'3pi/4_hi': vals}, open('ckpt_Mk_3pi4_hi.json','w'))
with open('progress.log','a') as f:
    f.write("step_60: M_k(3pi/4) k=0..3 recomputed at 300 dps quadrature (280 digits stored)\n")

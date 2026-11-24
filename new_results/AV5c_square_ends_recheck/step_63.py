import sympy as sp, mpmath, json, pickle, time
from mpmath import mp, mpf, pi, nstr
mp.dps = 300
u, g = sp.symbols('u g', positive=True)
ca = sp.symbols('ca')
fks = pickle.load(open('ckpt_corner_integrand.pkl','rb'))['fks']
cval = mp.cos(3*mp.pi/4)   # cos(alpha), alpha = 3pi/4
t0 = time.time()
vals = []
for k, fk in enumerate(fks):
    fkl = sp.lambdify((g, ca), fk, 'mpmath')
    fun = lambda x: fkl(x, cval)
    v = mpmath.quad(fun, [mpf(0), mpf('0.25'), mpf('0.5'), mpf('0.75'), mpf(1)])
    vals.append(nstr(v, 280))
    print(f"k={k}: {nstr(v,60)}   ({time.time()-t0:.0f}s)")
json.dump({'3pi/4_hi': vals}, open('ckpt_Mk_3pi4_hi.json','w'))
# cross-check k=0,1
c8 = mp.cos(3*mp.pi/8)
M0f = -(1+c8)*mp.log(1+c8)-(1-c8)*mp.log(1-c8)
print("k=0 vs elementary:", nstr(mpf(vals[0])-M0f, 6))
with open('progress.log','a') as f:
    f.write("step_61: M_k(3pi/4) k=0..3 at 300 dps with CORRECT ca=cos(3pi/4); k=0 matches elementary formula\n")

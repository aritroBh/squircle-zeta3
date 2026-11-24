import sympy as sp, mpmath, pickle, json
from mpmath import mp, mpf, sqrt, pi, zeta
mp.dps = 70
s = sp.symbols('s', positive=True)
f = {k: sp.sympify(v) for k,v in pickle.load(open('ckpt_fk_exprs.pkl','rb')).items()}

# numerical integration of each f_k over [0, inf)
mpmath.mp.dps = 70
def integ(expr):
    func = sp.lambdify(s, expr, 'mpmath')
    return mpmath.quad(func, [0, mpmath.inf])

F = {}
for k in range(1,9):
    F[k] = integ(f[k])
    print(f"F_{k} = {mp.nstr(F[k],50)}")

print()
for k in range(1,9):
    print(f"pred c{k+1} = 4sqrt2 F{k} = {mp.nstr(4*sqrt(2)*F[k],50)}")
pickle.dump({k: mp.nstr(v,60) for k,v in F.items()}, open('ckpt_Fk.pkl','wb'))

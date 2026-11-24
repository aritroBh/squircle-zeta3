import sympy as sp, mpmath, pickle
mpmath.mp.dps = 70
from mpmath import mpf, sqrt, zeta, nstr
s = sp.symbols('s', positive=True)
f = {k: sp.sympify(v) for k,v in pickle.load(open('ckpt_fk_exprs.pkl','rb')).items()}

def integ(expr, S=mpf(50)):
    func = sp.lambdify(s, expr, 'mpmath')
    return mpmath.quad(func, [0, 1, 4, 10, 20, S])

F = {}
for k in range(1,9):
    F[k] = integ(f[k])
    print(f"F_{k} = {nstr(F[k],45)}")
print()
preds = {}
for k in range(1,9):
    preds[k+1] = 4*sqrt(2)*F[k]
    print(f"pred c{k+1} = {nstr(preds[k+1],45)}")
pickle.dump({k: nstr(v,60) for k,v in F.items()}, open('ckpt_Fk.pkl','wb'))
json_ok = {str(k): nstr(v,60) for k,v in preds.items()}
import json; json.dump(json_ok, open('ckpt_pred_c_squircle.json','w'))
# eta-part checks
print()
print("check F1 = zeta(2)/4:", nstr(zeta(2)/4,45))
print("check F2 = zeta(3)/4:", nstr(zeta(3)/4,45))

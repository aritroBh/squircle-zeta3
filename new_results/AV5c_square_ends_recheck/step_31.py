import sympy as sp, json
u, g = sp.symbols('u g', positive=True)
ca = sp.symbols('ca', positive=True)  # cos alpha
lg = sp.log(g); l1g = sp.log(1+g)
S = 1 + g**2*sp.exp(-2*u*lg) - 2*ca*g*sp.exp(-u*lg)
A = sp.exp(-(1+u)*l1g) * sp.sqrt(S)
F = (A - 1) * sp.exp(u*lg) / g
ser = sp.series(F, u, 0, 4)
fks = [sp.expand(ser.coeff(u, k)) for k in range(4)]
for k, fk in enumerate(fks):
    print(k, "terms:", sp.count_ops(fk))
sp.save = True
import pickle
pickle.dump({'fks': fks}, open('ckpt_corner_integrand.pkl','wb'))
print("saved integrand coefficients")
with open('progress.log','a') as f:
    f.write("step_31: corner integrand [A-1]e^{u log g}/g expanded in u to order 3 (sympy), saved\n")

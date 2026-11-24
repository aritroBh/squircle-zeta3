import sympy as sp, mpmath, pickle
mpmath.mp.dps = 70
s = sp.symbols('s', positive=True)
f = {k: sp.sympify(v) for k,v in pickle.load(open('ckpt_fk_exprs.pkl','rb')).items()}
f1 = f[1]
print("f1 =", f1)
func = sp.lambdify(s, f1, 'mpmath')
for x in ['0.0','0.5','1.0','5.0','20.0','50.0']:
    print(x, mpmath.mp.nstr(func(mpmath.mpf(x)),20))
# expected: g1 = log(1+e^{-2s})
for x in ['0.5','1.0','5.0','20.0','50.0']:
    xx = mpmath.mpf(x)
    print('exact', x, mpmath.mp.nstr(mpmath.log(1+mpmath.exp(-2*xx)),20))

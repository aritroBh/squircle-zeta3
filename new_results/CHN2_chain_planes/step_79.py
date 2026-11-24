
import pickle, sympy as sp, numpy as np
forms10 = pickle.load(open('forms10.pkl','rb'))
a1,a2,a3,a4,a5 = sp.symbols('a1 a2 a3 a4 a5')
AV=[a1,a2,a3,a4,a5]
def to_sym(p):
    terms=[sp.Rational(c.numerator,c.denominator)*sp.Mul(*[AV[i]**e[i] for i in range(5) if e[i]]) for e,c in p.items()]
    return sp.Add(*terms) if terms else sp.Integer(0)
Sig = a1+a2+a3+a4+a5
C4 = sp.expand(sp.div(to_sym(forms10[8][3]), a2*a3*a4*Sig, *AV)[0])*10
Q7p = sp.expand(sp.div(to_sym(forms10[10][6]), a2*a3*a4*Sig, *AV)[0])*224
C4l = sp.lambdify(AV, C4, 'numpy')
Q7l = sp.lambdify(AV, Q7p, 'numpy')
rng = np.random.default_rng(7)
best=[]
C4p = sp.Poly(C4.subs(a3,1), a5)
for trial in range(20000):
    v1,v2,v4 = rng.uniform(0.02, 8, 3)
    coeffs = [float(c) for c in sp.Poly(C4p.subs({a1:v1,a2:v2,a4:v4}), a5).all_coeffs()]
    for r in np.roots(coeffs):
        if abs(r.imag)<1e-9 and r.real>1e-4:
            q = Q7l(v1,v2,1.0,v4,r.real)
            sc = (v1+v2+1+v4+r.real)**6
            best.append((abs(q)/sc, v1,v2,v4,r.real, q))
best.sort()
print("smallest |Q7'|/scale on {C4=0}:")
for b in best[:15]:
    print(f"  relerr={b[0]:.3e}  a=({b[1]:.4f},{b[2]:.4f},1,{b[3]:.4f},{b[4]:.4f})  Q7={b[5]:.3g}")

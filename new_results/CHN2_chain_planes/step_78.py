
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

# scan: normalize a3=1; grid over (a1,a2,a4) positive; solve C4=0 for a5>0 via polynomial roots
C4p = sp.Poly(C4.subs(a3,1), a5)
rng = np.random.default_rng(1)
found=[]
signs=[]
for trial in range(4000):
    v1,v2,v4 = rng.uniform(0.05, 6, 3)
    coeffs = [float(c) for c in sp.Poly(C4p.subs({a1:v1,a2:v2,a4:v4}), a5).all_coeffs()]
    rts = np.roots(coeffs)
    for r in rts:
        if abs(r.imag) < 1e-9 and r.real > 1e-6:
            a5v = r.real
            q = Q7l(v1,v2,1.0,v4,a5v)
            scale = max(1.0, abs(q))
            found.append((v1,v2,1.0,v4,a5v,q))
            signs.append(np.sign(q))
found = np.array(found)
print("positive points on {C4=0}:", len(found))
qs = found[:,5]
print("Q7' sign distribution: +:", np.sum(qs>0), " -:", np.sum(qs<0))
print("min |Q7'|/|a|^6:", np.min(np.abs(qs)/np.sum(found[:,:5],axis=1)**6))
print("Q7' range:", qs.min(), qs.max())

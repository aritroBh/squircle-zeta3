import sympy as sp, numpy as np, time
a1,a2,a4=sp.symbols('a1 a2 a4')
G=sp.sympify([l.split(':',1)[1].strip() for l in open('R_factors.txt') if 'deg 22' in l][0])
A_expr=sp.sympify(open('subres1_A.txt').read()); B_expr=sp.sympify(open('subres1_B.txt').read())
Gp=sp.Poly(G,a4)
deg_a4=Gp.degree()
print("G degree in a4:", deg_a4)
# coefficient functions c_k(a1,a2) for a4^k, k=0..deg
cf=[sp.lambdify((a1,a2), sp.expand(Gp.coeff_monomial(a4**k)), 'numpy') for k in range(deg_a4+1)]
# leading coeff (a4^16) as polynomial in a1,a2:
print("leading coeff a4^16 (in a1,a2):", sp.expand(Gp.coeff_monomial(a4**deg_a4)).free_symbols, "terms", len(sp.Add.make_args(sp.expand(Gp.coeff_monomial(a4**deg_a4)))))
A_l=sp.lambdify((a1,a2,a4), A_expr, 'numpy'); B_l=sp.lambdify((a1,a2,a4), B_expr, 'numpy')
# quick test: roots at (1,1)
coeffs=np.array([cf[k](1.0,1.0) for k in range(deg_a4,-1,-1)], dtype=float)
print("coeffs finite:", np.all(np.isfinite(coeffs)))
rts=np.roots(coeffs)
pos=[r.real for r in rts if abs(r.imag)<1e-6 and r.real>0]
print("numpy positive roots at (1,1):", sorted(pos))
# a5 and Q11 at these
Q11=pickle.load(open('Q11_poly.pkl','rb')) if False else None
import pickle
Q11=pickle.load(open('Q11_poly.pkl','rb'))
def q11_np(a1v,a2v,a3v,a4v,a5v):
    A_=[a1v,a2v,a3v,a4v,a5v]; tot=0.0
    for e,c in Q11.items():
        tot+=float(c)*np.prod([A_[i]**e[i] for i in range(5)])
    return tot
for a4r in sorted(pos):
    a5v=-B_l(1.0,1.0,a4r)/A_l(1.0,1.0,a4r)
    if a5v>0:
        print(f"  a4={a4r:.6f} a5={a5v:.6f} Q11={q11_np(1,1,1,a4r,a5v):.6g}")

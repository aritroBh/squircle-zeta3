import sympy as sp, numpy as np
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
# coefficient functions in a5 (C4 deg4, Q7 deg6), with a3=1, depending on (a1,a2,a4)
C4p=sp.Poly(C4.subs(a3,1), a5); Q7p=sp.Poly(Q7.subs(a3,1), a5)
print("C4 deg in a5:", C4p.degree(), " Q7 deg in a5:", Q7p.degree())
c4cf=[sp.lambdify((a1,a2,a4), sp.expand(C4p.coeff_monomial(a5**k)), 'numpy') for k in range(C4p.degree()+1)]
q7cf=[sp.lambdify((a1,a2,a4), sp.expand(Q7p.coeff_monomial(a5**k)), 'numpy') for k in range(Q7p.degree()+1)]
def common_pos_a5(a1v,a2v,a4v, tol=1e-6):
    """positive a5 that are roots of BOTH C4 and Q7' at (a1,a2,a4)."""
    cc=np.array([c4cf[k](a1v,a2v,a4v) for k in range(C4p.degree(),-1,-1)],dtype=float)
    if not np.all(np.isfinite(cc)): return []
    rc4=np.roots(cc)
    out=[]
    qv=np.array([q7cf[k](a1v,a2v,a4v) for k in range(Q7p.degree(),-1,-1)],dtype=float)
    for r in rc4:
        if abs(r.imag)<1e-7 and r.real>tol:
            # check Q7' at this a5
            val=np.polyval(qv, r.real)
            nrm=max(abs(qv)) 
            if abs(val)/nrm<1e-4:
                out.append(r.real)
    return out
# test at (1,1) for each positive a4 root
import numpy as np
G=sp.sympify([l.split(':',1)[1].strip() for l in open('R_factors.txt') if 'deg 22' in l][0])
Gp=sp.Poly(G,a4); deg_a4=Gp.degree()
gcf=[sp.lambdify((a1,a2), sp.expand(Gp.coeff_monomial(a4**k)), 'numpy') for k in range(deg_a4+1)]
coeffs=np.array([gcf[k](1.0,1.0) for k in range(deg_a4,-1,-1)],dtype=float)
for r in np.roots(coeffs):
    if abs(r.imag)<1e-6 and r.real>0:
        a5s=common_pos_a5(1.0,1.0,r.real)
        print(f"a4={r.real:.6f}: common positive a5 = {a5s}")

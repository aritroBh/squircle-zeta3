
import sympy as sp, numpy as np, pickle
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
G=sp.sympify([l.split(':',1)[1].strip() for l in open('R_factors.txt') if 'deg 22' in l][0])
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Gp=sp.Poly(G,a4); deg_a4=Gp.degree()
gcf=[sp.lambdify((a1,a2), sp.expand(Gp.coeff_monomial(a4**k)), 'numpy') for k in range(deg_a4+1)]
C4p=sp.Poly(C4.subs(a3,1), a5); Q7p=sp.Poly(Q7.subs(a3,1), a5)
d4,d6=C4p.degree(),Q7p.degree()
c4cf=[sp.lambdify((a1,a2,a4), sp.expand(C4p.coeff_monomial(a5**k)), 'numpy') for k in range(d4+1)]
q7cf=[sp.lambdify((a1,a2,a4), sp.expand(Q7p.coeff_monomial(a5**k)), 'numpy') for k in range(d6+1)]
Q11=pickle.load(open('Q11_poly.pkl','rb'))
QE=np.array(list(Q11.keys())); QC=np.array([float(c) for c in Q11.values()])
def q11_np(pt):
    pt=np.asarray(pt,dtype=float)
    pw=np.prod(np.power(pt[:,None], QE.T), axis=0)
    return float(np.dot(QC, pw))
def roots_pos(coeffs):
    coeffs=np.asarray(coeffs,dtype=float)
    if not np.all(np.isfinite(coeffs)): return []
    nz=np.argmax(np.abs(coeffs)>0); coeffs=coeffs[nz:]
    if len(coeffs)<2: return []
    m=np.max(np.abs(coeffs))
    if m==0: return []
    coeffs=coeffs/m
    try: r=np.roots(coeffs)
    except Exception: return []
    return [x.real for x in r if abs(x.imag)<1e-6 and x.real>1e-9]
def common_pos_a5(a1v,a2v,a4v):
    cc=np.array([c4cf[k](a1v,a2v,a4v) for k in range(d4,-1,-1)])
    out=[]
    for a5r in roots_pos(cc):
        qv=np.array([q7cf[k](a1v,a2v,a4v) for k in range(d6,-1,-1)])
        if abs(np.polyval(qv/np.max(np.abs(qv)), a5r))<1e-5: out.append(a5r)
    return out
def v10_points(a1v,a2v):
    cc=np.array([gcf[k](a1v,a2v) for k in range(deg_a4,-1,-1)])
    out=[]
    for a4r in roots_pos(cc):
        for a5r in common_pos_a5(a1v,a2v,a4r):
            out.append((a4r,a5r))
    return out

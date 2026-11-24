import sys, os, pickle
sys.path.insert(0, '<project>/research/08_theory/ALPHA2_angle_leading_law')
os.chdir('<project>/research/08_theory/ALPHA2_angle_leading_law')
from fractions import Fraction as Fr
from math import factorial
import sympy as sp

def qpart(t): return Fr(*t[0])
def fdiff(seq, n):
    for _ in range(n): seq = [seq[i+1]-seq[i] for i in range(len(seq)-1)]
    return seq
def extract_beta18(name):
    data = pickle.load(open(f'kern18_{name}.pkl','rb'))
    K = data['kernels']; r = Fr(*data['r']); c = Fr(*data['c'])
    tab = {}
    for d in range(1, 9):
        keys = [(i, j) for i in range(d, (d-1)//2, -1) for j in [d-i]]
        for key in keys:
            seqs = {}
            for k in range(d, 19):
                P = K[k][key]
                for rho, cf in P['real'].items(): seqs.setdefault(rho, []).append(qpart(cf))
            for rho, s in seqs.items():
                if len(s) < d + 3: continue
                top = fdiff(list(s[1:]), d-1)
                assert all(x == top[0] for x in top), (name, d, key, rho)
                beta = top[0]/factorial(d-1)
                if beta != 0: tab[(d, key, Fr(*rho))] = beta
    return c, r, tab

bs = pickle.load(open('beta_squircle.pkl','rb'))
sq = {}
for (d,ij,p),v in bs.items(): sq[(d,ij,{'0':0,'1':1,'-1':-1,'-r':2}[p])] = v[0]
Fd = pickle.load(open('Fd_data.pkl','rb'))
for name in ['c12','c35','c13','c45','c23','c57','c14','c67']:
    c, r, tab = extract_beta18(name)
    for d in [8]:
        ratios = set()
        for (dd,ij,p),v in tab.items():
            if dd!=d: continue
            marker = 2 if p == -r else int(p)
            if (d,ij,marker) in sq and sq[(d,ij,marker)]!=0: ratios.add(v/sq[(d,ij,marker)])
        if len(ratios)==1: Fd.setdefault(8,{})[c] = ratios.pop()
        else: print("issue", name, d, ratios)
print("d=8:", {str(c): Fd[8][c] for c in sorted(Fd[8])})
pickle.dump(Fd, open('Fd_data.pkl','wb'))
# fit F_8(t)
t = sp.symbols('t')
pts = [(sp.Rational(1), sp.Rational(1))]
for c, v in Fd[8].items():
    x = sp.Rational(c.numerator, c.denominator)
    pts.append(((1-x**2)/x**2, sp.Rational(v.numerator, v.denominator)))
for deg in range(3, 6):
    coeffs = sp.symbols(f'a0:{deg+1}')
    poly = sum(coeffs[i]*t**i for i in range(deg+1))
    sol = sp.solve([sp.Eq(poly.subs(t,tp), yp) for tp,yp in pts[:deg+1]], coeffs, dict=True)
    if not sol: continue
    p = sp.expand(poly.subs(sol[0]))
    if all(sp.simplify(p.subs(t,tp)-yp)==0 for tp,yp in pts):
        F = sp.together(p); N = sp.denom(F)
        print(f"d=8: F_8(t) = {p} | f_8 = {sp.expand(F*N)}, f_8(1)={N} ({len(pts)} pts, deg {deg})")
        fits = pickle.load(open('Fd_fits.pkl','rb')); fits[8] = p
        pickle.dump(fits, open('Fd_fits.pkl','wb'))
        break

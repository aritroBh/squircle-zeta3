import sys, os, pickle
sys.path.insert(0, '<project>/research/08_theory/ALPHA2_angle_leading_law')
os.chdir('<project>/research/08_theory/ALPHA2_angle_leading_law')
from assemble import beta_table
from fractions import Fraction as Fr
import sympy as sp
# squircle table has d<=8; need sq entries for d=7
bs = pickle.load(open('beta_squircle.pkl','rb'))
sq = {}
for (d,ij,p),v in bs.items(): sq[(d,ij,{'0':0,'1':1,'-1':-1,'-r':2}[p])] = v[0]
names = ['c12','c35','c13','c45','c23','c57','c14','c67']
Fd = pickle.load(open('Fd_data.pkl','rb'))
for name in names:
    c, r, tab = beta_table(name)
    for d in [7]:
        ratios = set()
        for (dd,ij,p),v in tab.items():
            if dd!=d: continue
            marker = 2 if p == -r else int(p)
            if (d,ij,marker) in sq and sq[(d,ij,marker)]!=0: ratios.add(v/sq[(d,ij,marker)])
        if len(ratios)==1: Fd.setdefault(d,{})[c] = ratios.pop()
        else: print("issue", name, d, ratios)
print("d=7 values:", {str(c): Fd[7][c] for c in sorted(Fd.get(7,{}))})
pickle.dump(Fd, open('Fd_data.pkl','wb'))
# fit d=7
t = sp.symbols('t')
pts = [(sp.Rational(1), sp.Rational(1))]
for c, v in Fd[7].items():
    x = sp.Rational(c.numerator, c.denominator)
    pts.append(((1-x**2)/x**2, sp.Rational(v.numerator, v.denominator)))
for deg in range(2, 6):
    coeffs = sp.symbols(f'a0:{deg+1}')
    poly = sum(coeffs[i]*t**i for i in range(deg+1))
    sol = sp.solve([sp.Eq(poly.subs(t,tp), yp) for tp,yp in pts[:deg+1]], coeffs, dict=True)
    if not sol: continue
    p = sp.expand(poly.subs(sol[0]))
    if all(sp.simplify(p.subs(t,tp)-yp)==0 for tp,yp in pts):
        F = sp.together(p); N = sp.denom(F)
        print(f"d=7: F_7(t) = {p}  | f_7 = {sp.expand(F*N)}, f_7(1)={N}  ({len(pts)} pts, deg {deg})")
        fits = pickle.load(open('Fd_fits.pkl','rb')); fits[7] = p
        pickle.dump(fits, open('Fd_fits.pkl','wb'))
        break

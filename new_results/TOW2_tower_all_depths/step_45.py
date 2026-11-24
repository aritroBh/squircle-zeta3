import sys, os, pickle
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
d, m = sp.symbols('d m')
i = sp.I
fits = pickle.load(open('e_fits.pkl','rb'))
def show(delta, s):
    r = fits[(delta,s)]
    mons = r['mons']; sol = r['sol']
    poly = sum(sp.Rational(re)*sp.nsimplify(im)*0 + sp.Rational(re)*d**ab[0]*m**ab[1] + sp.I*sp.Rational(im)*d**ab[0]*m**ab[1] for ab,(re,im) in zip(mons,sol))
    den = sp.prod(2*d-2*j-1 for j in range(1, r['b']+1)) if r['b'] else 1
    print(f'e({delta},{s})/e00 = ({sp.expand(poly)})/({den})   [b={r["b"]}]')
    return poly, den
p41, d41 = show(4,1)
# my hand-derived form
g = 4*d**3-36*d**2+113*d-116
mine = g/4 + sp.I*(2*m-d)/2
print('e(4,1) matches hand form:', sp.expand(p41 - mine) == 0)
p73, d73 = show(7,3)
zig = 5*sp.I*(d-5)/sp.prod(2*d-2*j-1 for j in (1,2,3))
print('e(7,3): fitted numerator:')
print(sp.expand(p73))
print('zigzag form:', sp.expand(zig*d73))
print('equal after *den?:', sp.expand(p73 - zig*d73) == 0)

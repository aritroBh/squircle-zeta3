import sys, pickle, os
work = '<project>/research/08_theory/TOW3R_review/'
sys.path.insert(0, work+'tow3_copy'); os.chdir(work+'tow3_copy')
import sympy as sp
import tow3lib as T3
d, m = T3.d, T3.m; i = sp.I
ef = pickle.load(open('e_fits.pkl','rb'))
def fit_expr(rec):
    b = rec['b']; mons = rec['mons']; sol = rec['sol']
    num = sum((re_+im_*i)*d**di*m**mi for (di,mi),(re_,im_) in zip(mons,sol))
    return sp.cancel(num / T3.Pi(b))
for key in [(4,1),(5,1),(6,2)]:
    fe = fit_expr(ef[key]); lib = T3._EF[key]
    t = sp.together(fe-lib)
    print(key, repr(t)[:120], '| simplify==0:', sp.simplify(t)==0, '| cancel==0:', sp.cancel(fe-lib)==0)

work = '<project>/research/08_theory/TOW3R_review/'
import sys, pickle
sys.path.insert(0, work+'tow3_copy'); import os; os.chdir(work+'tow3_copy')
import sympy as sp
import tow3lib as T3
d, m = T3.d, T3.m; i = sp.I
ef = pickle.load(open('e_fits.pkl','rb'))
def fit_expr(rec):
    b = rec['b']; mons = rec['mons']; sol = rec['sol']
    num = sum((re_+im_*i)*d**di*m**mi for (di,mi),(re_,im_) in zip(mons,sol))
    return sp.cancel(num / T3.Pi(b))
f41 = fit_expr(ef[(4,1)])
print('fit (4,1):', f41)
print('lib (4,1):', sp.cancel(T3.e41))
print('diff:', sp.expand(f41 - T3.e41))

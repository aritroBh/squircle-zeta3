import sys, os, pickle
work = '<project>/research/08_theory/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
import tow2lib as T
efits = pickle.load(open('e_fits.pkl','rb'))
dsym, msym = T.d, T.m
def fit_to_sym(f):
    num = sp.Integer(0)
    for (ii,jj),(re,im) in zip(f['mons'], f['sol']):
        num += (sp.Rational(re) + sp.Rational(im)*sp.I) * dsym**ii * msym**jj
    return num / T.Pi(f['b'])
e41_fit = sp.cancel(fit_to_sym(efits[(4,1)]))
e41_known = (dsym**3 - 9*dsym**2 + sp.Rational(113,4)*dsym - sp.I*dsym/2 + sp.I*msym - 29)/((2*dsym-3)*(2*dsym-5))
print('fit (4,1) == known:', sp.simplify(e41_fit - e41_known) == 0)
print(efits[(4,1)]['sol'])

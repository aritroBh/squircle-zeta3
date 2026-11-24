import sys, os, time, pickle
from fractions import Fraction as Fr
from math import comb
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
for key in [(5,1),(6,2)]:
    f = efits[key]
    expr = fit_to_sym(f)
    print(key, 'b=', f['b'], 'deg=', f['deg'])
    print('  ', sp.cancel(expr))

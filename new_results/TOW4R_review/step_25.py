
import os, sys, pickle
import sympy as sp
work = "<project>/research/08_theory/TOW4R_review"
d, m = sp.symbols('d m'); i = sp.I
ef = pickle.load(open(work + "/tow4_copy/e_fits.pkl", 'rb'))
v = ef[(7,3)]
print(type(v), v if not isinstance(v, dict) else {kk: (str(vv)[:120]) for kk, vv in v.items()})

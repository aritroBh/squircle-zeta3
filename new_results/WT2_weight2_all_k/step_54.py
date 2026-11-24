import sys, importlib, pickle, os
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
import fastred; importlib.reload(fastred)
from fastred import *

record = pickle.load(open('record2.pkl','rb'))
v = sp.Symbol('v')
S2 = sp.sqrt(2)
def real_form(Gm):
    return sp.together(pf_to_sympy(Gm))
g20_rec = real_form(record[(5,(2,0))]['Gm'])
g20_guess  = (1+S2)/4 - S2/(2*(v**2+1))
d = sp.together(g20_rec - g20_guess)
print("g20 diff:", [sp.N(d.subs(v, x), 20) for x in (sp.Rational(7,10), sp.Rational(13,10), 3)])
g11_2_rec = real_form(record[(2,(1,1))]['Gm'])
g11_2_guess = S2/(2*(v**2+1)) - (S2+1)/4
d2 = sp.together(g11_2_rec - g11_2_guess)
print("g11_2 diff:", [sp.N(d2.subs(v, x), 20) for x in (sp.Rational(7,10), sp.Rational(13,10), 3)])
print("g11_2_rec =", sp.expand(g11_2_rec))
print("guess expanded:", sp.expand(g11_2_guess))

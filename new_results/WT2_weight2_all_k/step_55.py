import sys, importlib, pickle, os
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
import fastred; importlib.reload(fastred)
from fastred import *
from fastred import v   # the module-level symbol (positive=True)

record = pickle.load(open('record2.pkl','rb'))
S2 = sp.sqrt(2)
def real_form(Gm):
    return sp.together(pf_to_sympy(Gm))
g20_rec = real_form(record[(5,(2,0))]['Gm'])
g11_2_rec = real_form(record[(2,(1,1))]['Gm'])
g11_3_rec = real_form(record[(5,(1,1))]['Gm'])
g20_guess  = (1+S2)/4 - S2/(2*(v**2+1))
g11_2_guess = S2/(2*(v**2+1)) - (S2+1)/4
g11_3_guess = -S2*(v-1)/(2*(v**2+1)) - S2/4
print("g20 match:", sp.simplify(sp.together(g20_rec-g20_guess))==0)
print("g11_2 match:", sp.simplify(sp.together(g11_2_rec-g11_2_guess))==0)
print("g11_3 match:", sp.simplify(sp.together(g11_3_rec-g11_3_guess))==0)
print("g20(v)  =", sp.factor(g20_rec))
print("g11_2(v)=", sp.factor(g11_2_rec))
print("g11_3(v)=", sp.factor(g11_3_rec))

import sys, importlib, pickle, os
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
import fastred; importlib.reload(fastred)
from fastred import *

record = pickle.load(open('record2.pkl','rb'))
v, eps = sp.symbols('v eps')
S2 = sp.sqrt(2); a = 1+S2; b = S2-1; I = sp.I

# Extract exact g20, g11_2, g11_3 from record (as real rational functions of v)
def real_form(Gm):
    return sp.together(pf_to_sympy(Gm))

g20_rec = real_form(record[(5,(2,0))]['Gm'])   # stable since k=2
g11_2_rec = real_form(record[(2,(1,1))]['Gm'])
g11_3_rec = real_form(record[(5,(1,1))]['Gm']) # stable since k=3
print("g20(v)  =", sp.simplify(g20_rec))
print("g11_2(v)=", sp.simplify(g11_2_rec))
print("g11_3(v)=", sp.simplify(g11_3_rec))

# my hand-derived guesses:
g20_guess  = (1+S2)/4 - S2/(2*(v**2+1))
g11_2_guess = S2/(2*(v**2+1)) - (S2+1)/4
g11_3_guess = -S2*(v-1)/(2*(v**2+1)) - S2/4
print("g20 match:", sp.simplify(g20_rec-g20_guess)==0)
print("g11_2 match:", sp.simplify(g11_2_rec-g11_2_guess)==0)
print("g11_3 match:", sp.simplify(g11_3_rec-g11_3_guess)==0)

# stabilization check: Gm_{20} same for all k in 2..12; Gm_{11} same for k in 3..12
stab20 = all(sp.simplify(real_form(record[(k,(2,0))]['Gm'])-g20_rec)==0 for k in range(2,13))
stab11 = all(sp.simplify(real_form(record[(k,(1,1))]['Gm'])-g11_3_rec)==0 for k in range(3,13))
print("Gm20 stable k=2..12:", stab20, " Gm11 stable k=3..12:", stab11)

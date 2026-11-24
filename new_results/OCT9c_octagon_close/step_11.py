import pickle, sys
import sympy as sp
d='../OCT9b_octagon_debug/'
def L(name):
    return pickle.load(open(d+name,'rb'))

brackets = L('brackets.pkl')
print(type(brackets), list(brackets.keys()) if hasattr(brackets,'keys') else len(brackets))

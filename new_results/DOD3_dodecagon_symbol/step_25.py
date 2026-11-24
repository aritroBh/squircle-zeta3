import sys, os, pickle
sys.path.insert(0, os.getcwd())
base = os.path.dirname(os.getcwd())
d2 = os.path.join(base,'DOD2_dodecagon_finish')
R = pickle.load(open(os.path.join(d2,'R256_sym.pkl'),'rb'))
import sympy as sp
print("R(2,5,6) as expression:")
print(R)

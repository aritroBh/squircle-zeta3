import pickle
from fractions import Fraction as Fr
base2 = "<project>/research/08_theory/Q1S2_vertex_close/"
symX = pickle.load(open(base2+"symX_V3p.pkl",'rb'))['T0']
VS = pickle.load(open(base2+"V_symbols.pkl",'rb'))['V3p'][0]
mine = pickle.load(open('q3_symbols.pkl','rb'))['S3_V3'][0]
print("mine (%d):"%len(mine))
for k,c in sorted(mine.items()): print("  ",k,c)
print("symX T0 (%d):"%len(symX))
for k,c in sorted(symX.items()): print("  ",k,c)
print("VS V3p[0] (%d):"%len(VS))
for k,c in sorted(VS.items()): print("  ",k,c)

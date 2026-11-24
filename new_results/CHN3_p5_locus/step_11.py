import pickle, sys
from fractions import Fraction as F
forms = pickle.load(open('forms10.pkl','rb'))
print("weights:", sorted(forms.keys()))
for W in sorted(forms.keys()):
    print(f"W={W}: {len(forms[W])} basis coords; nonzero:", 
          [b for b in range(len(forms[W])) if forms[W][b]])
# W=8 basis: [z8, z3z5, z2z3^2, z(3,5)] -> index 3 is z(3,5)
P8 = forms[8][3]
print("\nP8 terms:", len(P8))
# W=10 basis: [z10,z3z7,z5z5,z2z3z5,z4z3z3,z2z(3,5),z(3,7)] -> index 5,6
P10_235 = forms[10][5]
P10_37 = forms[10][6]
print("P10 z2z(3,5) terms:", len(P10_235), " z(3,7) terms:", len(P10_37))

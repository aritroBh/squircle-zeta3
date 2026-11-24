
import sys, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import math

reds = pickle.load(open('I_reductions_w2_7.pkl','rb'))
e_rm = pickle.load(open('e_rm_ext.pkl','rb'))
def E(r, m): return Fr(e_rm.get(f"{r},{m}", '0'))
bnames7 = ['log2^7', 'z2*log2^5', 'z2^2*log2^3', 'z2^3*log2', 'z3*log2^4', 'z3*z2*log2^2', 'z3*z2^2', 'z3^2*log2',
     'z5*log2^2', 'z5*z2', 'Li4*log2^3', 'Li4*z2*log2', 'Li4*z3', 'Li5*log2^2', 'Li5*z2', 'Li6*log2',
     'X6*log2', 'z7', 'Li7', 'X7a', 'X7b']

# for D~_{7,0}: contributions of each non-zeta element from each I
n, m = 7, 0
print(f"Assembly of D~_{{7,0}} = 2 sum_j e_{{0,6-j}} I_{{6-j,j}}/j! :")
print(f"{'element':>14} | " + " | ".join(f"I_{{{6-j},{j}}}" for j in range(1, 7)) + " | sum")
e_w = [2*E(0, 6-j)/math.factorial(j) for j in range(1, 7)]
print("e-weights:", [str(x) for x in e_w])
for i, nm in enumerate(bnames7):
    row = [reds[(6-j, j)][i]*e_w[j-1] for j in range(1, 7)]
    tot = sum(row)
    if tot or any(row):
        print(f"{nm:>14} | " + " | ".join(str(x) if x else "0" for x in row) + f" | {tot}")

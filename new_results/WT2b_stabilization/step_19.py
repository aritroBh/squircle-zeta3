import pickle
from fractions import Fraction as Fr
record = pickle.load(open('record_W2b.pkl','rb'))
# check order-1 +-i coefficient of every W_ij vanishes
bad = []
for (k,key), rec in record.items():
    W = rec['W']
    if W[1] and not W[1][0].is_zero(): bad.append((k,key,'+i',W[1][0]))
    if W[2] and not W[2][0].is_zero(): bad.append((k,key,'-i',W[2][0]))
print("nonzero order-1 +-i terms in any W_ij (k<=12):", bad if bad else "NONE -- exactness holds")

# check also: simple part of W_ij has no +-i at all; and print the surviving simple real-pole kernels at each key for k=3
rec = record[(3,(2,1))]
print("example W_21(k=3) real simple poles:", {str(r.sympy()):str(c.sympy()) for r,c in rec['W'][0].items()})

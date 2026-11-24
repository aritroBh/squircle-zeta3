
import pickle
from fractions import Fraction
PSI = pickle.load(open('/tmp/psi.pkl','rb'))['PSI']
# N_k = (-log2)^(k+1)/(k+1)! + sum_{m=0}^{k-1} INT[m], INT[m] = sum d_{ijk} (log2)^k J_{j,i}
def show_N(k):
    print(f"--- N_{k}:  (-log2)^{k+1}/{k+1}!  + ")
    for m in range(0, k):
        if not PSI[m]: continue
        for (i,j,kk),v in sorted(PSI[m].items()):
            print(f"      {v} * (log2)^{kk} J_{{{j},{i}}}")
show_N(2)
show_N(3)
show_N(4)

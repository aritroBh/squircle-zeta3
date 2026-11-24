import sys, os, time, pickle
from fractions import Fraction as Fr
from math import comb
work = '<project>/research/08_theory/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
import tow3lib as T3
H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
def u(i,j,k,delta):
    return tuple(H[(i,j)][delta][k]) if (i,j) in H and delta in H[(i,j)] and k in H[(i,j)][delta] else (Fr(0),)*4
# inspect key (2,2): u values at k=5,6,7 delta=1
for k in [4,5,6,7]:
    print('k',k, [u(2,2,k,dl) for dl in range(0,5)])
# Is (2,2) at k=5 perhaps off-law in the data (k=n+1 exception)? Check polynomiality of the data itself:
for dl in range(0,5):
    seq = [u(2,2,k,dl) for k in range(5,13)]
    # finite differences at order floor(dl/2)+1 should vanish
    P = dl//2
    cur = seq
    for _ in range(P+1):
        cur = [tuple(a-b for a,b in zip(cur[t+1],cur[t])) for t in range(len(cur)-1)]
    print('delta',dl,'(2,2) data poly check (should be all zero if law holds from k=5):', all(all(c==0 for c in cell) for cell in cur))

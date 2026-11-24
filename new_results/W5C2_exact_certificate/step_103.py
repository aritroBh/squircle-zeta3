import os, pickle, numpy as np
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
d = pickle.load(open('work/step01_pp.pkl','rb'))
S, C = d['S'], d['C']
feed = pickle.load(open('work/feed_pp.pkl','rb'))
print('|S|', len(S), '|C|', len(C), 'unique S', len(set(S)), 'unique C', len(set(C)))
# check triangularity of echelon pivots: rebuild piv rows? we didn't save them.
# Instead: rank of R_{S,C} mod p1 via LeanEchelon
import sys; sys.path.insert(0,'.')
from w5c2_lib import LeanEchelon, P1
p = P1
Cset = sorted(C); cpos = {c:k for k,c in enumerate(Cset)}
rows = []
for si in S:
    r = feed[si][2]
    cc = [cpos[c] for c,a in r if c in cpos]
    vv = [int(a)%p for c,a in r if c in cpos]
    rows.append((np.array(cc,np.int32), np.array(vv,np.int64)))
E = LeanEchelon(len(C), p)
import time
t0=time.monotonic()
E.feed(rows, logevery=4000, t0=t0)
print('rank of R_{S,C} mod p1 =', E.r, 'of', len(S), f'({time.monotonic()-t0:.0f}s)')
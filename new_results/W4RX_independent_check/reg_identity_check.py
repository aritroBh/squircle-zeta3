"""W4RX: exact check of the shuffle identity  w = sum_{k=0}^{m} (1)^{sh k}/k! sh reg(w[k:])
(G notation: letter 0 is the point 1; m = number of leading letters 0) for every word of
length 1..4 without trailing letter 8 (the point 0).  Used in the proof of the regularised
distribution relation."""
import sys
sys.dont_write_bytecode = True
from math import factorial
from itertools import product
from w4rx_rel import *
bad = 0; n = 0
for L in range(1, 5):
    for w in product(range(9), repeat=L):
        if w[-1] == ZERO: continue
        m = 0
        while m < len(w) and w[m] == 0: m += 1
        tot = dict(reg(w)); Ck = {(): F(1)}
        for k in range(1, m + 1):
            Ck = vmul(Ck, {(0,): F(1)})
            tot = vadd((tot, 1), (vmul(Ck, reg(w[k:])), F(1, factorial(k))))
        n += 1
        bad += tot != {w: 1}
print('words', n, 'failures', bad)

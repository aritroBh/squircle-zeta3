import sys, os, pickle
sys.path.insert(0,'.')
import mpmath as mp
from flint import fmpz_mat
mp.mp.dps = 300
# replicate with a few atoms, including L32's from kpool
kpool = pickle.load(open('kpool_1550.pkl','rb'))
z3 = mp.zeta(3); z5 = mp.zeta(5); z2 = mp.zeta(2)
vals_d = {'z5': z5, 'z2z3': z2*z3, 'L32pp': mp.mpf(kpool['L32pp']['K']), 'L32qq': mp.mpf(kpool['L32qq']['K'])}
INV = list(vals_d)
n = len(INV); r5 = mp.sqrt(5)
D = 250
C = mp.mpf(10)**int(D*0.92)
rows = []
for i in range(n):
    row = [0]*(n+2); row[i] = 1
    row[n] = int(mp.nint(C*vals_d[INV[i]])); row[n+1] = int(mp.nint(C*r5*vals_d[INV[i]]))
    rows.append(row)
M = fmpz_mat(rows).lll()
for i in range(n):
    v = [int(M[i,j]) for j in range(n)]
    tail = [int(M[i,j]) for j in range(n, n+2)]
    res1 = abs(sum(mp.mpf(v[j])*vals_d[INV[j]] for j in range(n)))
    res5 = abs(sum(mp.mpf(v[j])*r5*vals_d[INV[j]] for j in range(n)))
    print("row", i, "h", max(abs(x) for x in v), "resid", mp.nstr(res1,2), mp.nstr(res5,2), "tail", tail)

import sys, os, pickle
sys.path.insert(0,'.')
import mpmath as mp
from flint import fmpz_mat
mp.mp.dps = 200
kpool = pickle.load(open('kpool_1550.pkl','rb'))
z2 = mp.zeta(2); z3 = mp.zeta(3); z5 = mp.zeta(5); r5 = mp.sqrt(5)
phi = (1+r5)/2; lf = mp.log(phi); l2 = mp.log(2)
vals_d = {'z2l2^3': z2*l2**3, 'z3lf^2': z3*lf**2, 'z5': z5, 'L32pp': mp.mpf(kpool['L32pp']['K']), 'z2^2l2': z2**2*l2}
INV = list(vals_d); n = len(INV)
for cexp in (150,):
    C = mp.mpf(10)**cexp
    rows = []
    for i in range(n):
        row = [0]*(n+2); row[i] = 1
        row[n] = int(mp.nint(C*vals_d[INV[i]])); row[n+1] = int(mp.nint(C*r5*vals_d[INV[i]]))
        rows.append(row)
    print("input row 0:", rows[0][-2:], " row 4:", rows[4][-2:])
    M = fmpz_mat(rows).lll()
    for i in range(n):
        v = [int(M[i,j]) for j in range(n)]
        tail = [int(M[i,j]) for j in range(n, n+2)]
        res1 = abs(sum(mp.mpf(v[j])*vals_d[INV[j]] for j in range(n)))
        # verify lattice membership: tail[0] should equal sum v[j]*rows[j][n]
        t0 = sum(v[j]*rows[j][n] for j in range(n))
        t1 = sum(v[j]*rows[j][n+1] for j in range(n))
        print("row", i, "h %.2e" % max(abs(x) for x in v), "resid1", mp.nstr(res1,2), "tail", tail, "expected", (t0,t1), "match:", tail==[t0,t1])

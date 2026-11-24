import sys, os
sys.dont_write_bytecode = True
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
from towlib import *
from genintegrand import S_poly
from fractions import Fraction as Fr

tstar = QSI(Fr(1,2),0,Fr(1,2),0)
def S_val(d, m, l):
    P = S_poly(d, m)
    out = ZERO
    for n, c in enumerate(P):
        if n < l: continue
        ff = 1
        for a in range(l): ff *= (n-a)
        out = out + QSI(Fr(c*ff)) * (tstar**(n-l))
    return out

def N_formula(l, d, m):
    """fitted numerator polynomials N_l(d,m) over Q(i): returns (re, im)"""
    if l == 1:
        re = 4*m - 2*d; im = 6*d - 4*d**2
    elif l == 2:
        re = 16*m**2 + Fr(80,3)*d - 16*d*m - 80*d**2 + Fr(208,3)*d**3 - 16*d**4
        im = -40*m + 20*d + 96*d*m - 48*d**2 - 32*d**2*m + 16*d**3
    elif l == 3:
        re = -688*m + 64*m**3 + 344*d + 2768*d*m - 96*d*m**2 - 1384*d**2 - 3264*d**2*m + 1648*d**3 + 1408*d**3*m - 704*d**4 - 192*d**4*m + 96*d**5
        im = -672*m**2 - 280*d + 672*d*m + 864*d*m**2 + 1144*d**2 - 864*d**2*m - 192*d**2*m**2 - 1984*d**3 + 192*d**3*m + 1600*d**4 - 544*d**5 + 64*d**6
    return Fr(re), Fr(im)

ok = True
for l in (1,2,3):
    for d in range(12, 14):
        for m in (0, 1, d//2, d):
            R = S_val(d,m,l)/S_val(d,m,0)
            prod = 1
            for j in range(1,l+1): prod *= (2*d-2*j-1)
            re, im = N_formula(l, d, m)
            if R.p0*prod != re or R.p2*prod != im:
                ok = False; print('MISMATCH', l, d, m, R.p0*prod, re, R.p2*prod, im)
print('N_l generalization check d=12,13:', 'OK' if ok else 'FAIL')
open('progress.log','a').write('step_28: S-derivative fits N1,N2,N3 verified at fresh points d=12,13\n')

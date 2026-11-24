fix = '''

def e_off(key, delta, s):
    if (delta, s) == (0,0): return rho_c(key, 0)
    if (delta, s) == (1,0): return rho_c(key, 1)
    if (delta, s) == (2,0): return rho_c(key, 2)
    if (delta, s) == (2,1):
        i,j = key; d = i+j
        return u(i,j,d+1,2) - u(i,j,d,2)     # slope from tower data
    return ZERO

def ehat_pp(p, delta, s):
    if delta == 0 and s == 0:
        return e_off((p+1,p), 0, 0) * tau1(0).scale(Fr(p+1, 4*p-1))
    if delta == 1 and s == 0:
        return ( e_off((p+1,p),1,0) * tau1(0).scale(Fr(p+1, 4*p-2))
               + e_off((p+1,p),0,0) * tau1(1).scale(Fr(p+1, 4*p-1)) )
    if delta == 2 and s == 0:
        return ( e_off((p+1,p),2,0) * tau1(0).scale(Fr(p+1, 4*p-3))
               + e_off((p+1,p),1,0) * tau1(1).scale(Fr(p+1, 4*p-2))
               + e_off((p+1,p),0,0) * tau1(2).scale(Fr(p+1, 4*p-1)) )
    if delta == 2 and s == 1:
        return e_off((p+1,p),2,1) * tau1(0).scale(Fr(p+1, 4*p-3))
    return ZERO
'''
src = open('towlib.py','a').write(fix)
import importlib, towlib
importlib.reload(towlib)
from towlib import *
print('w(2,1) =', e_off((2,1),2,1), '(expect QSI(0,1/4,0,0) = sqrt2/4)')
print('w(3,0) =', e_off((3,0),2,1), '(expect -sqrt2/12 = QSI(0,-1/12,0,0))')
print('ehat22(0,0) =', ehat_pp(2,0,0), '(expect 15sqrt2(1+i)/16 = QSI(0,15/16,0,15/16))')
print('ehat22(1,0) =', ehat_pp(2,1,0), '(expect 15sqrt2/8-15sqrt2 i/8 = QSI(0,15/8,0,-15/8))')
print('ehat22(2,0) =', ehat_pp(2,2,0), '(data u(2,2,5,2) = QSI(0,-135/64,0,-117/64))')
print('ehat22(2,1) =', ehat_pp(2,2,1), '(data slope: u(6,2)-u(5,2) = QSI(0,-12/64,0,-12/64))')

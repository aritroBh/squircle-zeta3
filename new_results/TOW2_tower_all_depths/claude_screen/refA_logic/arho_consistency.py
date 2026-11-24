import sys, time
sys.path.insert(0,'.')
import sympy as sp
import tow2lib as T
d,m,i = T.d,T.m,sp.I
t0=time.time()
# A. Bell-assembled arho (used by NA proofs) vs pickle arho_simple (used by AD/E3 proofs)
for delta in (1,2,3,4):
    diff = sp.cancel(sp.together(T.arho(delta) - T.arho_simple(delta)))
    print('arho vs arho_simple, delta=',delta,'->',diff, round(time.time()-t0,1))
# TOW1's rho(2)/rho(0) = 1/2[(3d-3)/2 - (i/2) r_p - (d-2)^2], r_p = 2(d-2m) - 2 d i - 4(d-1)(d-2m)/(2d-3)
r = 2*(d-2*m) - 2*d*i - 4*(d-1)*(d-2*m)/(2*d-3)
rho2_TOW1 = sp.Rational(1,2)*((3*d-3)/sp.Integer(2) - i*r/2 - (d-2)**2)
print('TOW1 rho(2) vs arho_simple(2):', sp.cancel(sp.together(rho2_TOW1 - T.arho_simple(2))))
print('TOW1 rho(1) = -i(d-2):', sp.cancel(T.arho_simple(1) - (-i*(d-2))))
# proof-doc a3, a4 text vs pickle
a3 = (4*i*d**4 - 36*i*d**3 + (6+125*i)*d**2 - 12*d*m + (-18-198*i)*d + 36*m + 117*i)/(24*d-36)
a4 = (16*d**6 - 240*d**5 + (1528-48*i)*d**4 + 96*i*d**3*m + (-5296+432*i)*d**3 - 864*i*d**2*m + (10525-1308*i)*d**2 + (48+2616*i)*d*m + (-11321+1320*i)*d - 48*m**2 - 2640*i*m + 5040)/(384*d**2-1536*d+1440)
print('proof-doc a3 vs pickle:', sp.cancel(sp.together(a3-T.arho_simple(3))))
print('proof-doc a4 vs pickle:', sp.cancel(sp.together(a4-T.arho_simple(4))))
# E. does arho(delta<=4) depend on N[3], N[4]? (sl_3, sl_4 enter only L_r with r>=6)
import inspect
print('L_r uses sl_l for l <= r//2; for r<=4 -> l<=2. sl keys:', list(T.sl.keys()))
# check symbolic: arho(4) free of N[3],N[4] by recomputing with garbage sl[3], sl[4]
old3, old4 = T.sl[3], T.sl[4]
T.sl[3] = sp.Symbol('Z3'); T.sl[4] = sp.Symbol('Z4'); T._L.clear(); T._a.clear()
a4g = T.arho(4); a3g = T.arho(3)
print('arho(3),arho(4) free of garbage sl3/sl4:', not (a4g.has(T.sl[3]) or a4g.has(T.sl[4]) or a3g.has(T.sl[3]) or a3g.has(T.sl[4])))

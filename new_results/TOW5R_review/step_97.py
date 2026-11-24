import sys
wd = '<project>/research/08_theory/TOW5R_review'
sys.path.insert(0, wd+'/tow5_copy')
import sympy as sp
import tow2lib as T2, tow5lib as T5
d, m = T2.d, T2.m
i = sp.I
Pi4 = T5.Pi(4); Pi5 = T5.Pi(5)
f93 = ( -20*i*d**4 + 420*i*d**3 + (-30 - 3457*i)*d**2 + 60*d*m + (180 + 12768*i)*d - 360*m - 16956*i ) / (12*Pi4)
f103 = ( -80*d**6 + 2640*d**5 + (-37496 + 240*i)*d**4 - 480*i*d**3*m + (290096 - 4320*i)*d**3
         + 8640*i*d**2*m + (-1271513 + 26508*i)*d**2 + (-240 - 53016*i)*d*m + (2937427 - 53532*i)*d
         + 240*m**2 + 107064*i*m - 2733468 ) / (96*Pi5)
print('F(9,3) display match:', sp.cancel(sp.together(T5.forced(9,3) - f93)) == 0)
print('F(10,3) display match:', sp.cancel(sp.together(T5.forced(10,3) - f103)) == 0)
# F(9,1), F(10,1), F(10,2) are "[deg ... terms]/(den Pi)" claims: check denominators only
n91, d91 = sp.fraction(sp.cancel(T5.forced(9,1)))
n101, d101 = sp.fraction(sp.cancel(T5.forced(10,1)))
n102, d102 = sp.fraction(sp.cancel(T5.forced(10,2)))
print('F(9,1): den =', sp.factor(d91), '; num deg_d', sp.degree(n91,d), 'terms', len(sp.Poly(n91,d,m).terms()))
print('F(10,1): den =', sp.factor(d101), '; num deg_d', sp.degree(n101,d), 'terms', len(sp.Poly(n101,d,m).terms()))
print('F(10,2): den =', sp.factor(d102), '; num deg_d', sp.degree(n102,d), 'terms', len(sp.Poly(n102,d,m).terms()))

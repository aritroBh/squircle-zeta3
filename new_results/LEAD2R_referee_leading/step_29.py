import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 60
from lead2r_lib import *

A1 = zigzag_ints(25)
A2 = zigzag_recurrence(25)
print("zigzag two ways agree:", A1 == A2)
print("A_0..10:", A1[:11])
ok = all(eps_m(m) == [1,-1,-1,1][m%4] for m in range(20))
print("eps formula == printed period-4:", ok)
# Lemma 6.1 first identity: delta_d = eps_{d-1} A_d/d!
from math import factorial
okd = all(delta_rat(d) == Fr(eps_m(d-1)*A1[d], factorial(d)) for d in range(1, 22))
print("delta_d = eps_{d-1} A_d/d! for d=1..21:", okd)
# second identity: eps_m A_{m+1}/m! = [z^m](1-sinh z)/cosh^2 z
import sympy as sp
z = sp.symbols('z')
ser = sp.series((1-sp.sinh(z))/sp.cosh(z)**2, z, 0, 22).removeO()
ok2 = all(sp.Rational(ser.coeff(z,m)) == sp.Rational(eps_m(m)*A1[m+1], factorial(m)) for m in range(0,21))
print("eps_m A_{m+1}/m! = [z^m](1-sinh)/cosh^2, m=0..20:", ok2)
# E coeffs two ways
E1 = e_coeffs_logE(20); E2 = e_coeffs_taylor(20)
worst = max(abs(E1[n]-E2[n]) for n in range(21))
print("E coeffs logE vs taylor worst diff:", mp.nstr(worst,3))
print("e_0..5:", [mp.nstr(E1[n],8) for n in range(6)])
with open('progress.log','a') as f:
    f.write("step_03: sanity of own lib: zigzag 2 ways agree; delta=eps A/d! d<=21; [z^m](1-sinh)/cosh^2 identity m<=20; E coeffs 2 ways worst "+mp.nstr(worst,3)+"\n")

"""v_basis_change.py -- basis dependence of the zeta(3) coefficient at 2pi/3 (verifier).
(1) the two classical identities that move between bases A, B, C (100 digits);
(2) the three W3 forms printed in m2_2pi3_pslq.out, typed by hand, + W2 = z2/2 - Li2(1/4)/2, vs solver's M2 (318 dig).
"""
import mpmath as mp
mp.mp.dps = 330
L3 = lambda x: mp.polylog(3, x); L2 = lambda x: mp.polylog(2, x)
l2, l3, pi2, z3 = mp.log(2), mp.log(3), mp.pi**2, mp.zeta(3)
x = mp.mpf(1)/3
id1 = L3(x) + L3(1 - x) + L3(1 - 1/x) - (z3 + mp.log(x)**3/6 + pi2/6*mp.log(x) - mp.log(x)**2*mp.log(1 - x)/2)
id2 = L3(-2) - L3(-mp.mpf(1)/2) - (-pi2/6*l2 - l2**3/6)
id3 = L3(mp.mpf(1)/2) + L3(-mp.mpf(1)/2) - L3(mp.mpf(1)/4)/4
id4 = L3(mp.mpf(1)/2) - (7*z3/8 - pi2*l2/12 + l2**3/6)
print('Landen x=1/3:', mp.nstr(id1, 3), ' inversion -2:', mp.nstr(id2, 3), ' duplication 1/2:', mp.nstr(id3, 3), ' Li3(1/2):', mp.nstr(id4, 3))
M2 = [mp.mpf(l.split('=')[1]) for l in open('values_2pi_3.txt') if l.startswith('M2 =')][0]
W2 = pi2/12 - L2(mp.mpf(1)/4)/2
common = lambda: (-mp.mpf(7)/2*L2(x)*l2 + mp.mpf(15)/4*L2(x)*l3)
A = mp.mpf(15)/2*L3(x) + mp.mpf(17)/2*L3(-mp.mpf(1)/2) + mp.mpf(3)/16*z3 - mp.mpf(5)/12*pi2*l2 + pi2*l3/4 - l2**3/12 \
    + mp.mpf(5)/4*l2**2*l3 - mp.mpf(7)/4*l2*l3**2 + mp.mpf(5)/8*l3**3 + common()
B = -mp.mpf(15)/2*L3(mp.mpf(2)/3) + L3(-mp.mpf(1)/2) + mp.mpf(123)/16*z3 + mp.mpf(5)/6*pi2*l2 - pi2*l3 + mp.mpf(7)/6*l2**3 \
    + mp.mpf(5)/4*l2**2*l3 - mp.mpf(11)/2*l2*l3**2 + mp.mpf(25)/8*l3**3 + common()
C = mp.mpf(15)/2*L3(x) + mp.mpf(17)/8*L3(mp.mpf(1)/4) - mp.mpf(29)/4*z3 + mp.mpf(7)/24*pi2*l2 + pi2*l3/4 - mp.mpf(3)/2*l2**3 \
    + mp.mpf(5)/4*l2**2*l3 - mp.mpf(7)/4*l2*l3**2 + mp.mpf(5)/8*l3**3 + common()
for nm, v in (('A (z3: 3/16)', A), ('B (z3: 123/16)', B), ('C (z3: -29/4)', C)):
    print(f'basis {nm}: |W2 + W3 - M2(2pi/3)| =', mp.nstr(abs(W2 + v - M2), 3))

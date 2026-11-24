import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from fractions import Fraction as F
from qsqrt2 import QS2 as Q
import mpmath as mp
mp.mp.dps = 100

rel1500 = {k: tuple(F(x) for x in v) for k, v in pickle.load(open("rel1500_frac.pkl","rb")).items()}
R2 = Q(0, 1)
def Q2(a, b): return Q(a, b)
# zeta(5) coefficient in P5/u
z5c = Q(rel1500['z5'][0], rel1500['z5'][1])
# zeta(2)zeta(3) explicit + combined with pi2*z3 = 6 zeta(2)zeta(3)
z23a = Q(rel1500['z2*z3'][0], rel1500['z2*z3'][1])
z23b = Q(rel1500['z3*pi2'][0], rel1500['z3*pi2'][1])
z23_comb = z23a + 6*z23b
print("coeff of zeta(5) in P^(5)(1)/u:        (%s)" % z5c)
print("coeff of zeta(2)zeta(3) explicit:      (%s)" % z23a)
print("coeff of pi2*zeta(3):                  (%s)" % z23b)
print("coeff of zeta(2)zeta(3) [combined]:    (%s)" % z23_comb)

u = mp.sqrt(2+mp.sqrt(2)); r2 = mp.sqrt(2)
def q2f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r2
print("\nin P^(5)(1) = u*(P^(5)(1)/u):")
print("  zeta(5)        coeff =", mp.nstr(u*q2f(z5c), 30))
print("  zeta(2)zeta(3) coeff =", mp.nstr(u*q2f(z23_comb), 30), "(combined)")

# local law m=4
s2 = (2-R2)/4; s4 = s2*s2; s6 = s4*s2
A5 = 65*s6 - 63*s4 + 17*s2 - 1
A23 = 45*s6 - 27*s4 + 6*s2 - 1
C5 = 2880*4*s4*A5; C23 = -960*4*s4*A23
c = u/2
print("\nlocal law (prop:order5, m=4):")
print("  C_4^(5)   = (%s)*c = %s" % (C5, mp.nstr(q2f(C5)*c, 30)))
print("  C_4^(2,3) = (%s)*c = %s" % (C23, mp.nstr(q2f(C23)*c, 30)))
print("\nnon-local parts (actual - local), basis-dependent:")
print("  zeta(5):        %s" % mp.nstr(u*q2f(z5c) - q2f(C5)*c, 25))
print("  zeta(2)zeta(3): %s" % mp.nstr(u*q2f(z23_comb) - q2f(C23)*c, 25))

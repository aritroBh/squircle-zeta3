import sys
sys.path.insert(0, '<project>/research/08_theory/TOW4_tower_depth8')
import sympy as sp
import tow2lib as T2
import mpmath as mp
mp.mp.dps = 50
i = sp.I
s2 = mp.sqrt(2); aq = s2+1; bq = s2-1; I = mp.mpc(0,1)
def tau1_mp(n):
    s = n+1
    return (-1)**s * ((I-bq)**(-(s+1)) + (I+aq)**(-(s+1)) - I**(-(s+1)))
def tau2_mp(n):
    s = n+1
    return (-1)**s * ((I-aq)**(-(s+1)) + (I+bq)**(-(s+1)) - I**(-(s+1)))
v1 = sp.Rational(-169,128) + sp.Rational(41,128)*i
v2 = sp.Rational(-169,128) - sp.Rational(41,128)*i
print("tau1(8) numeric:", mp.nstr(tau1_mp(8), 30))
print("cand1 numeric :", mp.nstr(mp.mpc(float(sp.re(v1)), float(sp.im(v1))), 30))
print("match1:", abs(tau1_mp(8) - (mp.mpf(-169)/128 + mp.mpc(0,1)*mp.mpf(41)/128)) < mp.mpf(10)**(-40))
print("match2:", abs(tau2_mp(8) - (mp.mpf(-169)/128 - mp.mpc(0,1)*mp.mpf(41)/128)) < mp.mpf(10)**(-40))
# exact check: definition equals candidate symbolically
d1 = sp.simplify(sp.expand(T2.tau1(8) - v1))
d2 = sp.simplify(sp.expand(T2.tau2(8) - v2))
print("symbolic tau1(8)==cand:", d1==0, " tau2(8)==cand:", d2==0)
# conjugation symmetry conj(tau1_8) == (-1)^8 tau2_8
print("conj symmetry:", sp.conjugate(v1).subs(sp.conjugate(sp.I),-sp.I) == v2 or (sp.re(v1)-sp.im(v1)*i)==v2)

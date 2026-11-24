import os, sys, time, pickle, re
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
import mpmath as mp
mp.mp.dps = 120
import gengine
from flint import fmpz_mat
from fractions import Fraction as F

# Rebuild the column set WITH values and symbols, get exact kernel, test values against 13-basis
import symengine as S
from symengine import Q, ZERO, ONE, P_LET, Q_LET, add, scale, shuffle, const_pow, li_symbol, gsym
R2 = Q(0,1); B = Q(-1,1)
ZLIST = [('b',B),('-b',-B),('b^2',B*B),('-b^2',-B*B),('1/2',Q(F(1,2))),('1/sqrt2',R2/2),
         ('sqrt2',R2),('-sqrt2',-R2),('2-sqrt2',2-R2),('p',P_LET),('q',Q_LET),('-1',-ONE),('2',Q(2))]
zmap = {'b':R2-1,'-b':-(R2-1),'b^2':(R2-1)**2,'-b^2':-(R2-1)**2,'1/2':mp.mpf('0.5'),'1/sqrt2':1/R2,
        'sqrt2':R2,'-sqrt2':-R2,'2-sqrt2':2-R2,'p':2+R2,'q':-(1+R2),'-1':mp.mpf(-1),'2':mp.mpf(2)}
X = mp.log(2); Y = mp.log(1+R2)
def li_val(n, name):
    z = zmap[name]
    return mp.re(mp.polylog(n, z)) if abs(z) > 1 else mp.polylog(n, z)

cols = []  # (label, symbol, value-func)
for name,z in ZLIST:
    cols.append((f"Li5({name})", li_symbol(5, z), lambda n=name: li_val(5,n)))
for name,z in ZLIST:
    s4 = li_symbol(4, z)
    cols.append((f"Li4({name})*X", shuffle(s4, const_pow(0,1)), lambda n=name: li_val(4,n)*X))
    cols.append((f"Li4({name})*Y", shuffle(s4, const_pow(1,1)), lambda n=name: li_val(4,n)*Y))
for name,z in ZLIST:
    s3 = li_symbol(3, z)
    cols.append((f"Li3({name})*X^2/2", shuffle(s3, const_pow(0,2)), lambda n=name: li_val(3,n)*X**2/2))
    cols.append((f"Li3({name})*X*Y", shuffle(s3, shuffle(const_pow(0,1),const_pow(1,1))), lambda n=name: li_val(3,n)*X*Y))
    cols.append((f"Li3({name})*Y^2/2", shuffle(s3, const_pow(1,2)), lambda n=name: li_val(3,n)*Y**2/2))
for name,z in ZLIST:
    s2 = li_symbol(2, z)
    cols.append((f"Li2({name})*X^3/6", shuffle(s2, const_pow(0,3)), lambda n=name: li_val(2,n)*X**3/6))
    cols.append((f"Li2({name})*X^2Y/2", shuffle(s2, shuffle(const_pow(0,2),const_pow(1,1))), lambda n=name: li_val(2,n)*X**2*Y/2))
    cols.append((f"Li2({name})*XY^2/2", shuffle(s2, shuffle(const_pow(0,1),const_pow(1,2))), lambda n=name: li_val(2,n)*X*Y**2/2))
    cols.append((f"Li2({name})*Y^3/6", shuffle(s2, const_pow(1,3)), lambda n=name: li_val(2,n)*Y**3/6))
for i in range(6):
    cols.append((f"X^{i}Y^{5-i}/({i}!{5-i}!)", shuffle(const_pow(0,i), const_pow(1,5-i)), lambda i=i: X**i*Y**(5-i)/(mp.factorial(i)*mp.factorial(5-i))))
for n1,z1 in ZLIST:
    for n2,z2 in ZLIST:
        cols.append((f"Li3({n1})*Li2({n2})", shuffle(li_symbol(3,z1), li_symbol(2,z2)), lambda a=n1,b2=n2: li_val(3,a)*li_val(2,b2)))
for n1,z1 in ZLIST:
    for n2,z2 in ZLIST:
        s22 = shuffle(li_symbol(2,z1), li_symbol(2,z2))
        cols.append((f"Li2({n1})*Li2({n2})*X", shuffle(s22, const_pow(0,1)), lambda a=n1,b2=n2: li_val(2,a)*li_val(2,b2)*X))
        cols.append((f"Li2({n1})*Li2({n2})*Y", shuffle(s22, const_pow(1,1)), lambda a=n1,b2=n2: li_val(2,a)*li_val(2,b2)*Y))
print("columns:", len(cols))
words = sorted(set().union(*[set(c[1]) for c in cols]))
M = [[c[1].get(w, ZERO) for c in cols] for w in words]
# exact kernel via sympy over QQ(sqrt2): use QS2 entries -> sympy
import sympy as sp
Ms = sp.Matrix([[sp.Rational(c.p.numerator, c.p.denominator) + sp.Rational(c.q.numerator, c.q.denominator)*sp.sqrt(2) for c in row] for row in M])
t0=time.time()
ns = Ms.nullspace()
print("kernel dim:", len(ns), "in %.1fs" % (time.time()-t0))

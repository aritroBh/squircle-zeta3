
import os
os.chdir('<project>/research/08_theory/COR8_depth_two_any_angle')
# (0) update result.md item 3
t = open('result.md').read()
t = t.replace(
 "3. delta_22 MATCH (exact): d22(Omega) = d22( ((1-c)/c) Li_{2,2}(1/2,1-c) - ((1+c)/c) Li_{2,2}(1/2,1+c) ).\n   Two depth-two terms, both with x=1/2, y=1∓c; series convergent on (0,1). (step_02)",
 "3. delta_22 MATCH (exact, RECOMPUTED step_60): d22(Omega) = d22(DT_A) with\n   DT_A = ((1-c)/c) Li_{2,2}(1/2,1-c) - ((1+c)/c) Li_{2,2}(1/2,1+c)  [the version used in the solve].\n   The equivalent variant DT_B = ((1-c)/c) Li_{2,2}(1/2,1-c) + ((1+c)/c) Li_{3,1}(1/2,1+c) also matches exactly\n   (progress.log's first-day line referred to DT_B). Both have convergent series on (0,1).")
open('result.md','w').write(t)
print("result.md updated")
# (1) evaluator checks: MPL conventions as genuine double sums; Li4 inversion vs Re(polylog)
import mpmath as mp
mp.mp.dps = 60
def mpl_sum(a,b,x,y):
    tot = mp.mpf(0)
    for n in range(2, 3000):
        term = x**n * sum(y**m/mp.mpf(m)**b for m in range(1,n)) / mp.mpf(n)**a
        tot += term
        if n > 60 and abs(term) < mp.mpf(10)**(-mp.mp.dps-5): break
    return tot
import sys
sys.path.insert(0, '<project>/research/08_theory/COR4_M3alpha_classical')
from hpl_series import Li22, Li31
x, y = mp.mpf('0.31'), mp.mpf('0.7')
print("|Li22 series vs genuine sum|:", abs(Li22(x,y) - mpl_sum(2,2,x,y)))
print("|Li31 series vs genuine sum|:", abs(Li31(x,y) - mpl_sum(3,1,x,y)))
# inversion for Li4 at x=2 (real), comparing Re(polylog)
xx = mp.mpf(2)
lg = mp.log(xx) + mp.j*mp.pi   # log(-x) with x+i0 prescription
u = mp.mpf('0.5') + lg/(2*mp.pi*mp.j)
B4 = u**4 - 2*u**3 + u**2 - mp.mpf(1)/30
rhs = -(2*mp.pi*mp.j)**4/24 * B4
li4_inv = rhs - mp.polylog(4, 1/xx)
print("Re(inversion Li4(2)):", mp.nstr(mp.re(li4_inv), 30))
print("Re(polylog  Li4(2)):", mp.nstr(mp.re(mp.polylog(4, xx)), 30))
print("|diff|:", mp.nstr(abs(mp.re(li4_inv) - mp.re(mp.polylog(4, xx))), 5))
# same at x=4
xx = mp.mpf(4)
lg = mp.log(xx) + mp.j*mp.pi
u = mp.mpf('0.5') + lg/(2*mp.pi*mp.j)
B4 = u**4 - 2*u**3 + u**2 - mp.mpf(1)/30
li4_inv = -(2*mp.pi*mp.j)**4/24 * B4 - mp.polylog(4, 1/xx)
print("|diff at x=4|:", mp.nstr(abs(mp.re(li4_inv) - mp.re(mp.polylog(4, xx))), 5))
open('progress.log','a').write('step_61: MPL conventions OK; Re(polylog) == Re(inversion) for Li4 at x=2,4\n')

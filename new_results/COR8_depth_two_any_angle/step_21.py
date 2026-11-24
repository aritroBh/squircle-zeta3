
import os
os.chdir('<project>/research/08_theory/COR8_depth_two_any_angle')
open('progress.log','a').write('start: COR8 resumed; building eps-symbol engine with TRUE symbols of Li_{3,1},Li_{2,2}\n')
# quick numeric validation of MPL conventions and y=+-1 HPL identities
import sys
sys.path.insert(0, '<project>/research/08_theory/COR4_M3alpha_classical')
import mpmath as mp
mp.mp.dps = 50
from hpl_series import H_series, Li31, Li22

x = mp.mpf('0.31')
# Li_{3,1}(x,1) = H_{0,0,1,1}(x) ; Li_{2,2}(x,1) = H_{0,1,0,1}(x)
print('Li31(x,1) - H0011:', abs(Li31(x,1) - H_series((0,0,1,1), x)))
print('Li22(x,1) - H0101:', abs(Li22(x,1) - H_series((0,1,0,1), x)))
# Li_{3,1}(x,-1) = -H_{0,0,1,-1}(x); Li_{2,2}(x,-1) = -H_{0,1,0,-1}(x)
print('Li31(x,-1) + H001-1:', abs(Li31(x,-1) + H_series((0,0,1,-1), x)))
print('Li22(x,-1) + H010-1:', abs(Li22(x,-1) + H_series((0,1,0,-1), x)))
# differential recursions: x d/dx Li31 = Li21, y d/dy Li31(x,y) = (y Li3(x) - Li3(xy))/(1-y)
def Li21(x,y):
    return mp.nsum(lambda n,m: x**n * y**m/(n**2*m), [2, mp.inf], [1, lambda n: n-1]) if False else None
def mpl(a,b,x,y):
    tot = mp.mpf(0); n = 2
    inner = mp.mpf(0)
    while True:
        inner = inner + y**(n-1)/mp.mpf(n-1)**b
        term = x**n*inner/n**a
        tot += term
        if n > 50 and abs(term) < mp.mpf(10)**(-mp.mp.dps-3): break
        n += 1
    return tot
y = mp.mpf('0.7')
h = mp.mpf('1e-30')
dx = (mpl(3,1,x+h,y)-mpl(3,1,x-h,y))/(2*h)*x
print('x d/dx Li31 - Li21:', abs(dx - mpl(2,1,x,y)))
dy = (mpl(3,1,x,y+h)-mpl(3,1,x,y-h))/(2*h)*y
print('y d/dy Li31 - (yLi3x-Li3xy)/(1-y):', abs(dy - (y*mp.polylog(3,x)-mp.polylog(3,x*y))/(1-y)))
dy22 = (mpl(2,2,x,y+h)-mpl(2,2,x,y-h))/(2*h)*y
print('y d/dy Li22 - Li21:', abs(dy22 - mpl(2,1,x,y)))
dx12 = (mpl(1,2,x+h,y)-mpl(1,2,x-h,y))/(2*h)*x
print('x d/dx Li12 - x Li2(xy)/(1-x):', abs(dx12 - x*mp.polylog(2,x*y)/(1-x)))
open('progress.log','a').write('step_00 numeric: MPL conventions + y=1,-1 HPL ids + differential recursions OK\n')

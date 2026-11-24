import sys; sys.dont_write_bytecode=True
from mpmath import mp, mpf, gamma, hyp3f2, quad, rf
mp.dps=40
def vol16(u1,u2,u3,u4):
    G=gamma(1+u2)*gamma(1+u1+u3)/gamma(1+u1+u2+u3)
    return G*hyp3f2(-u4,u3,1+u1+u3,1+u3,1+u1+u2+u3,1)
def vol16_int(u1,u2,u3,u4):
    # u2 u3 * int_{s,t>=0,s+t<=1} s^(u2-1) t^(u3-1) (1-s)^u1 (1-t)^u4
    f=lambda s: s**(u2-1)*(1-s)**u1 * quad(lambda t: t**(u3-1)*(1-t)**u4,[0,1-s])
    return u2*u3*quad(f,[0,1])
def G(x): return gamma(1+x)
def pi1(x,y):   # a=(x,x,y,y)
    return (G(x)**2*G(y)**2/(G(2*x)*G(2*y)) + G(x+y)**2/G(2*x+2*y))/2
def pi2(x,y):   # a=(x,y,x,y)
    return (G(2*x)*G(2*y)/G(2*x+2*y) + G(x)**2*G(y)**2/G(x+y)**2)/2
def lam_minus(u1,u2,u3,u4):   # F11 - F12 on the s23=0 slice (generic a) : Gamma(1+u1+u3)Gamma(1+u2+u4)/Gamma(1+sum)
    return G(u1+u3)*G(u2+u4)/G(u1+u2+u3+u4)
for (x,y) in [(mpf('0.37'),mpf('0.81')),(mpf('1.25'),mpf('0.4')),(mpf('2.3'),mpf('1.7')),(mpf('0.05'),mpf('3.1'))]:
    v1=vol16(x,x,y,y); v2=vol16(x,y,x,y)
    print('x=%s y=%s'%(x,y))
    print('  Pi1 (u1=u2=x,u3=u4=y): 3F2 form %s'%mp.nstr(v1,30),' closed form diff = %s'%mp.nstr(v1-pi1(x,y),5))
    print('  Pi2 (u1=u3=x,u2=u4=y): 3F2 form %s'%mp.nstr(v2,30),' closed form diff = %s'%mp.nstr(v2-pi2(x,y),5))
# double-integral check of the task's integral form at one Pi1 point and one Pi2 point (lower precision)
mp.dps=20
x,y=mpf('0.6'),mpf('1.3')
print('integral form Pi1:',mp.nstr(vol16_int(x,x,y,y),15),' closed',mp.nstr(pi1(x,y),15))
print('integral form Pi2:',mp.nstr(vol16_int(x,y,x,y),15),' closed',mp.nstr(pi2(x,y),15))
# negative control: generic point is NOT of that form
mp.dps=30
print('generic (0.3,0.5,0.7,0.9) vol16 =',mp.nstr(vol16(mpf('0.3'),mpf('0.5'),mpf('0.7'),mpf('0.9')),20))

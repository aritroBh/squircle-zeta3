
import mpmath as mp
mp.mp.dps = 50

def vol16(u1,u2,u3,u4):
    G = mp.gamma(1+u2)*mp.gamma(1+u1+u3)/mp.gamma(1+u1+u2+u3)
    h = mp.hyper([ -u4, u3, 1+u1+u3],[1+u3, 1+u1+u2+u3], 1)
    return G*h

def Bterm(u1,u2,u3,u4):
    return mp.gamma(1+u1)*mp.gamma(1+u2)*mp.gamma(1+u3)*mp.gamma(1+u4)/(mp.gamma(1+u1+u2)*mp.gamma(1+u3+u4))
def G0term(u1,u2,u3,u4):
    return mp.gamma(1+u1+u3)*mp.gamma(1+u2+u4)/mp.gamma(1+u1+u2+u3+u4)

# point on Pi2: (x,y)=(1,2) -> a=(1,2,1,2)
x,y=mp.mpf(1),mp.mpf(2)
# extract c2 numerically: f(u) = 1 + c2 u^2 + c3 u^3 + ... use Richardson on f at a few u
import numpy as np
def c2_numeric(scale=1):
    # use f(u) with u = h, h/2: c2 ~ (f(h)-1)/h^2 corrected
    hs=[mp.mpf(1)/1000, mp.mpf(1)/2000, mp.mpf(1)/4000]
    vals=[]
    for h in hs:
        f=vol16(x*h,y*h,x*h,y*h)
        vals.append((f-1)/h**2)
    # Richardson extrapolate (error ~ h)
    v1 = 2*vals[1]-vals[0]
    v2 = 2*vals[2]-vals[1]
    return 2*v2-v1  # second-order
c2n = c2_numeric()
z2 = mp.zeta(2)
print("numeric c2 at (1,2,1,2):", mp.nstr(c2n,30))
print("(B+G0)/2 prediction  -13/2*z2 :", mp.nstr(-mp.mpf(13)/2*z2,30))
print("engine prediction    -6*z2      :", mp.nstr(-6*z2,30))
# also direct numeric of (B+G0)/2 vs vol16 at one u
u=mp.mpf('0.013')
v=vol16(x*u,y*u,x*u,y*u)
bg=(Bterm(x*u,y*u,x*u,y*u)+G0term(x*u,y*u,x*u,y*u))/2
print("vol16 - (B+G0)/2 at u=0.013:", mp.nstr(v-bg,10))

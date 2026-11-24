
import mpmath as mp
mp.mp.dps = 80

def vol16(u1,u2,u3,u4):
    G = mp.gamma(1+u2)*mp.gamma(1+u1+u3)/mp.gamma(1+u1+u2+u3)
    h = mp.hyper([ -u4, u3, 1+u1+u3],[1+u3, 1+u1+u2+u3], 1)
    return G*h
def bg(u1,u2,u3,u4):
    B = mp.gamma(1+u1)*mp.gamma(1+u2)*mp.gamma(1+u3)*mp.gamma(1+u4)/(mp.gamma(1+u1+u2)*mp.gamma(1+u3+u4))
    G0= mp.gamma(1+u1+u3)*mp.gamma(1+u2+u4)/mp.gamma(1+u1+u2+u3+u4)
    return (B+G0)/2

x,y=mp.mpf(1),mp.mpf(2)
for h in [mp.mpf(1)/100, mp.mpf(1)/1000, mp.mpf(1)/10000]:
    fv=vol16(x*h,y*h,x*h,y*h); fb=bg(x*h,y*h,x*h,y*h)
    print(f"h={mp.nstr(h,6)}  (fv-1)/h^2={mp.nstr((fv-1)/h**2,25)}")
    print(f"          (fb-1)/h^2={mp.nstr((fb-1)/h**2,25)}  fv-fb={mp.nstr(fv-fb,5)}")
print("-13/2*zeta(2) =", mp.nstr(-mp.mpf(13)/2*mp.zeta(2),25))
print("-6*zeta(2)    =", mp.nstr(-6*mp.zeta(2),25))

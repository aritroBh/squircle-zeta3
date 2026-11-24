
import mpmath as mp
mp.mp.dps = 60

def vol16(u1,u2,u3,u4):
    G = mp.gamma(1+u2)*mp.gamma(1+u1+u3)/mp.gamma(1+u1+u2+u3)
    h = mp.hyper([ -u4, u3, 1+u1+u3],[1+u3, 1+u1+u2+u3], 1)
    return G*h
def Bterm(u1,u2,u3,u4):
    return mp.gamma(1+u1)*mp.gamma(1+u2)*mp.gamma(1+u3)*mp.gamma(1+u4)/(mp.gamma(1+u1+u2)*mp.gamma(1+u3+u4))
def G0term(u1,u2,u3,u4):
    return mp.gamma(1+u1+u3)*mp.gamma(1+u2+u4)/mp.gamma(1+u1+u2+u3+u4)

x,y=mp.mpf(1),mp.mpf(2); u=mp.mpf('0.013')
v=vol16(x*u,y*u,x*u,y*u)
bg=(Bterm(x*u,y*u,x*u,y*u)+G0term(x*u,y*u,x*u,y*u))/2
print("vol16:", mp.nstr(v,40))
print("bg   :", mp.nstr(bg,40))
print("diff :", mp.nstr(v-bg,10))
# predicted diff from c2 analysis: (x-y)^2/2*z2*u^2
print("predicted diff:", mp.nstr((x-y)**2/2*mp.zeta(2)*u**2, 10))
# Also test the planner's point (0.3,0.7,0.3,0.7)
u1,u2,u3,u4 = [mp.mpf(t) for t in ('0.3','0.7','0.3','0.7')]
v2=vol16(u1,u2,u3,u4); bg2=(Bterm(u1,u2,u3,u4)+G0term(u1,u2,u3,u4))/2
print("planner pt diff:", mp.nstr(v2-bg2,10))

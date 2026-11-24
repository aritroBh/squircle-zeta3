
import mpmath as mp
mp.mp.dps = 60
alpha=mp.pi/2; c=mp.cos(alpha)
L=lambda x: mp.log1p(mp.exp(-x))
def Tint(x,q):
    z=x*(q-1); ez=mp.exp(z)
    return mp.exp(-q*L(x))*ez*(mp.sqrt(1-2*c*ez+ez**2)-1)/(mp.exp(x)+1)
def Rint(x,q):
    return (mp.exp(-q*L(x))/(1+mp.exp(-x))-1)*mp.exp(-q*x)
for q in [mp.mpf(0), mp.mpf('0.3')]:
    T1=mp.quad(lambda x: Tint(x,q),[0,1,10,mp.inf])
    R1=mp.quad(lambda x: Rint(x,q),[0,1,10,mp.inf])
    # finer quadrature
    T2=mp.quad(lambda x: Tint(x,q),[0,mp.mpf('0.5'),2,20,mp.inf])
    R2=mp.quad(lambda x: Rint(x,q),[0,mp.mpf('0.5'),2,20,mp.inf])
    print("q=%s T1=%s R1=%s"%(mp.nstr(q,4),mp.nstr(T1,20),mp.nstr(R1,20)))
    print("      T2=%s R2=%s  T2-T1=%.1e R2-R1=%.1e"%(mp.nstr(T2,20),mp.nstr(R2,20),abs(T2-T1),abs(R2-R1)))
    f=lambda g: ((1+g)**(-1-q)*mp.sqrt(1+g**(2-2*q)-2*c*g**(1-q))-1)*g**(q-1)
    M=mp.quad(f,[0,1])
    print("   M=%s  T+R=%s  diff=%.2e"%(mp.nstr(M,20),mp.nstr(T2+R2,20),abs(M-T2-R2)))

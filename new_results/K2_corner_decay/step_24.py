
import mpmath as mp
mp.mp.dps = 60

def M_direct(q, alpha):
    c=mp.cos(alpha)
    f=lambda g: ((1+g)**(-1-q)*mp.sqrt(1+g**(2-2*q)-2*c*g**(1-q))-1)*g**(q-1)
    return mp.quad(f,[0,mp.mpf('0.01'),mp.mpf('0.5'),mp.mpf('0.99'),1])

def M_TR(q, alpha):
    c=mp.cos(alpha)
    L=lambda x: mp.log1p(mp.exp(-x))
    def Tintegrand(x):
        z=x*(q-1); ez=mp.exp(z)
        Bt=mp.exp(-z)*(mp.sqrt(1-2*c*ez+ez**2)-1)   # FIXED: e^{-z}
        return mp.exp(-q*L(x))*Bt/(mp.exp(x)+1)
    def Rintegrand(x):
        return (mp.exp(-q*L(x))/(1+mp.exp(-x))-1)*mp.exp(-q*x)
    return mp.quad(Tintegrand,[0,1,10,mp.inf]) + mp.quad(Rintegrand,[0,1,10,mp.inf])

for alpha in [mp.pi/2, mp.pi/3, 5*mp.pi/6]:
    for q in [mp.mpf('0.3'), mp.mpf('0.7'), mp.mpc('0.3','0.7')]:
        v1=M_direct(q,alpha); v2=M_TR(q,alpha)
        print("alpha=%.4f q=%s  diff=%.2e" % (float(alpha), mp.nstr(q,6), abs(v1-v2)))

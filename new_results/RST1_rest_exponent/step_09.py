"""RST1 step_09: corrected Gaussian+H F-lemma v2.
Split [0, rb], [rb, rt], [rt, inf), rb = rm/2, rt = 2 r*, rm = r* - c0/2.
left  <= (2^{3/2}/lam_L) e^{E~(rb)}        [lam_L = E~'(rb), direct eval]
mid   <= x_b^{-3/2} sqrt(2 pi/c2) e^{-sqrt(pi j q)+ a c0/2}
right <= x_t^{-3/2} e^{-sqrt(pi j q)}/a
F_H <= sqrt(pi/2) (sqrt2/q) (left+mid+right).
"""
import mpmath as mp
mp.mp.dps = 50
pi=mp.pi; omega=mp.expjpi(mp.mpf(1)/4); a=pi/(2*mp.sqrt(2)); sq=mp.sqrt(2)
def H(s): return mp.beta(s/2, mp.mpf(3)/2)/2
def FHexact(j,q):
    s=lambda r: mp.mpf('0.5')+r*omega
    rstar=mp.sqrt(2*j*q/pi)
    f=lambda r: mp.e**(-a*r)*abs(H(s(r)))*abs(s(r))**j/abs(s(r)+q)**(j+1)
    pts=[0,rstar/4,rstar/2,rstar,2*rstar,4*rstar,10*rstar,mp.inf]
    return mp.quad(f,pts)
def pieces(j,q):
    rs=mp.sqrt(2*j*q/pi); c0=1+q/sq; rm=rs-c0/2
    rb=rm/2; rt=2*rs
    c2=8*sq*j*q/(2*rt+c0)**3
    sp=mp.sqrt(pi*j*q)
    lamL=-a+2*sq*j*q/(2*rb+c0)**2
    Eb=-a*rb-sq*j*q/(2*rb+c0)
    xb=mp.mpf('0.5')+rb/sq; xt=mp.mpf('0.5')+rt/sq
    left=2**mp.mpf('1.5')/lamL*mp.e**Eb
    mid=xb**(-mp.mpf('1.5'))*mp.sqrt(2*pi/c2)*mp.e**(-sp+a*c0/2)
    right=xt**(-mp.mpf('1.5'))*mp.e**(-sp)/a
    return (left+mid+right)*sq/q*mp.sqrt(pi/2), left, mid, right
def bound(j,q):
    return pieces(j,q)[0]
worst=0
for i in range(0,13):
    for k in [99,150,200,300,400,1000,3000]:
        j=k-i; q=i+1
        Fe=FHexact(j,q); b=bound(j,q); ra=Fe/b
        worst=max(worst,ra)
print("worst ratio exact/bound:", mp.nstr(worst,6))
print("bound*sqrt(2 j q)*e^{sqrt(pi j q)} (want ~const):")
for k in [99,200,400,1000,3000,10000]:
    j,q=k,2
    b,le,mi,ri=pieces(j,q)
    print(f"  j={j:6d}: tot={mp.nstr(b*mp.sqrt(2*j*q)*mp.e**(mp.sqrt(pi*j*q)),8)}  (left,mid,right)/tot = {mp.nstr(le/b,4)},{mp.nstr(mi/b,4)},{mp.nstr(ri/b,4)}")

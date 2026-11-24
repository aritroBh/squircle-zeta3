"""RST1 step_08: numerical test of the Gaussian+H F-lemma (power j^(-1/2)).
F_H(j,q) = int_0^inf e^{-a r} |H(s)| |s|^j/|s+q|^{j+1} dr, s=1/2+r omega.
Bound (r_*=sqrt(2jq/pi) > c0=1+q/sqrt2):
  (sqrt2/q) e^{-sqrt(pi j q)} [ (r*/(2sqrt2))^(-3/2) sqrt(2pi/c2) e^{a c0/2}
     + 2^(3/2) (r*/2) e^{a c0/2 - sqrt(pi j q)/64} + (2r*/sqrt2)^(-3/2) /a ] * sqrt(pi/2)
c2 = 8 sqrt2 j q/(4 r*+c0)^3.
Test exact/bound ratio for the needed (j,q) = (k-i, i+1), i=0..12, k=99..400.
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
def bound(j,q):
    rs=mp.sqrt(2*j*q/pi); c0=1+q/sq
    assert rs>c0, (j,q)
    c2=8*sq*j*q/(4*rs+c0)**3
    sp=mp.sqrt(pi*j*q)
    mid=(rs/(2*sq))**(-mp.mpf('1.5'))*mp.sqrt(2*pi/c2)*mp.e**(a*c0/2)
    left=2**mp.mpf('1.5')*(rs/2)*mp.e**(a*c0/2-sp/64)
    right=(2*rs/sq)**(-mp.mpf('1.5'))/a
    return mp.sqrt(pi/2)*sq/q*mp.e**(-sp)*(mid+left+right)
worst=0; wq=None
for i in range(0,13):
    for k in [99,150,200,300,400]:
        j=k-i; q=i+1
        Fe=FHexact(j,q); b=bound(j,q); ra=Fe/b
        if ra>worst: worst=ra; wq=(i,k,j,q,Fe,b)
print("worst ratio exact/bound over tested grid:", mp.nstr(worst,6), "at (i,k,j,q)=",wq[0],wq[1],wq[2],wq[3])
print("sample values:")
for (j,q) in [(99,1),(99,2),(400,2),(400,13),(87,13)]:
    print(f"  j={j} q={q}: exact={mp.nstr(FHexact(j,q),6)} bound={mp.nstr(bound(j,q),6)}")
# also check the exponent prefactor scaling: bound * sqrt(j) e^{+sqrt(pi j q)}
print("bound*sqrt(2jq)*e^{sqrt(pi j q)} (should be ~bounded constant):")
for k in [99,200,400]:
    j,q=k,2
    print(f"  j={j}: {mp.nstr(bound(j,q)*mp.sqrt(2*j*q)*mp.e**(mp.sqrt(pi*j*q)),8)}")

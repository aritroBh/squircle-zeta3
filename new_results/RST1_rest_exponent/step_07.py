"""RST1 step_07: numerical test of the F-lemma (sharp-exponent modulus bound).
F(j,q) = int_0^inf e^{-a r} |s|^j/|s+q|^{j+1} dr,  s=1/2+r*omega, a=pi/(2 sqrt 2).
Bound: sqrt2/q * ( 2 sqrt(2 j q/pi) e^{a(1+q/sqrt2)/2} + 1/a ) e^{-sqrt(pi j q)}.
Test ratio exact/bound for q=2, j=10..400 and a few (j,q).
"""
import mpmath as mp
mp.mp.dps = 50
pi=mp.pi; omega=mp.expjpi(mp.mpf(1)/4); a=pi/(2*mp.sqrt(2)); sq=mp.sqrt(2)
def Fexact(j,q):
    s=lambda r: mp.mpf('0.5')+r*omega
    rstar=mp.sqrt(2*j*q/pi)
    f=lambda r: mp.e**(-a*r)*abs(s(r))**j/abs(s(r)+q)**(j+1)
    pts=[0,rstar/4,rstar/2,rstar,2*rstar,4*rstar,10*rstar,mp.inf]
    return mp.quad(f,pts)
def bound(j,q):
    return sq/q*(2*mp.sqrt(2*j*q/pi)*mp.e**(a*(1+q/sq)/2)+1/a)*mp.e**(-mp.sqrt(pi*j*q))
print(" q=2:   j      Fexact            bound            ratio")
for j in [10,50,100,200,400,1000]:
    Fe=Fexact(j,2); b=bound(j,2)
    print(f"  {j:5d}  {mp.nstr(Fe,6):>12s}  {mp.nstr(b,6):>12s}  {mp.nstr(Fe/b,6)}")
print(" other (j,q):")
for (j,q) in [(100,1),(100,3),(100,4),(50,9),(10,2),(3,20),(0,50),(1,101)]:
    Fe=Fexact(j,q); b=bound(j,q)
    print(f"  j={j:4d} q={q:3d}  {mp.nstr(Fe,6):>12s}  {mp.nstr(b,6):>12s}  ratio={mp.nstr(Fe/b,6)}")

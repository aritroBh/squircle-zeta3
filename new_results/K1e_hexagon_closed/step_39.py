
import mpmath as mp
mp.mp.dps=80
z2v=mp.pi**2/6
print("=== inversion (odd weight): Li3(-x) - Li3(-1/x) = -(1/6)ln^3 x - (pi^2/6) ln x ===")
for xv in ['2','3','1/2','7','1']:
    x=mp.mpf(xv)
    lhs=mp.polylog(3,-x)-mp.polylog(3,-1/x)
    rhs=-mp.log(x)**3/6 - z2v*mp.log(x)
    print(f"x={xv:>3}: |d|={mp.nstr(abs(lhs-rhs),3)}")

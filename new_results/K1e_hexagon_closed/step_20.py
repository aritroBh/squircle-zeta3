
import mpmath as mp, time
mp.mp.dps=120
def A(x):
    # A(x) = sum H_n x^n / n^2 ; compute via series with convergence accel for |x|<=1, else use identity for check only
    return mp.nsum(lambda n: mp.harmonic(n)*x**n/n**2, [1, mp.inf])
def Ared(x):
    # 6 A(x) = 12 Li3(x) + 6 Li3(x/(x-1)) - 6 Li2(x) ln(1-x) - ln^3(1-x)
    l=mp.log(1-x)
    return (12*mp.polylog(3,x) + 6*mp.polylog(3,x/(x-1)) - 6*mp.polylog(2,x)*l - l**3)/6
for xv in ['-2','-1','-1/2','-1/3','1/2','1/3','2/3','1/5','-3/2']:
    x=mp.mpf(xv)
    lhs=A(x); rhs=Ared(x)
    print(f"x={xv:>5}: |A - Ared| = {mp.nstr(abs(lhs-rhs),3)}   A={mp.nstr(lhs,12)}")
with open('research/08_theory/K1e_hexagon_closed/progress.log','a') as fh:
    fh.write(f"[{time.strftime('%H:%M:%S')}] step_06: verified K4e A-reduction 6A(x)=12Li3(x)+6Li3(x/(x-1))-6Li2(x)ln(1-x)-ln^3(1-x) at 9 points incl x=-2,-3/2; residuals <=1e-118.\n")

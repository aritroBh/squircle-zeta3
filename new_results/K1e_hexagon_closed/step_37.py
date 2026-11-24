
import mpmath as mp
mp.mp.dps=50
for xv in ['2/3','1/3']:
    x=mp.mpf(xv)
    print(f"x={xv}: Li3(x)={mp.nstr(mp.polylog(3,x),15)}, Li3(1-x)={mp.nstr(mp.polylog(3,1-x),15)}, Li3(1-1/x)={mp.nstr(mp.polylog(3,1-1/x),15)} (arg {mp.nstr(1-1/x,5)})")
    lhs=mp.polylog(3,x)+mp.polylog(3,1-x)+mp.polylog(3,1-1/x)
    # correct Lewin form candidates:
    c1 = mp.zeta(3) + mp.pi**2/6*mp.log(1-x) - mp.log(x)*mp.log(1-x)**2/2 + mp.log(1-x)**3/6
    c2 = mp.zeta(3) + mp.pi**2/6*mp.log(1-x) - mp.log(-x)*mp.log(1-x)**2/2 + mp.log(1-x)**3/6
    c3 = mp.zeta(3) - mp.pi**2/6*mp.log(1-x) - mp.log(x)*mp.log(1-x)**2/2 + mp.log(1-x)**3/6
    print(f"  lhs={mp.nstr(lhs,15)} c1={mp.nstr(c1,15)} |d1|={mp.nstr(abs(lhs-c1),3)}")
    print(f"     |d2(ln(-x))|={mp.nstr(abs(lhs-c2),3)}  |d3(-z2 ln)|={mp.nstr(abs(lhs-c3),3)}")

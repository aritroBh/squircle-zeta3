
import mpmath as mp
mp.mp.dps=110
z3v=mp.zeta(3); z2v=mp.pi**2/6
for xv in ['1/2','2/3','1/3','3/4']:
    x=mp.mpf(xv)
    t1=mp.polylog(3,x); t2=mp.polylog(3,1-x); t3=mp.polylog(3,1-1/x)
    lhs=t1+t2+t3
    r1=z3v; r2=z2v*mp.log(1-x); r3=-mp.log(x)*mp.log(1-x)**2/2; r4=mp.log(1-x)**3/6
    print(f"x={xv}: lhs={mp.nstr(lhs,15)}")
    print(f"   rhs pieces: z3={mp.nstr(r1,8)}, z2*ln(1-x)={mp.nstr(r2,8)}, -.5 ln x ln^2={mp.nstr(r3,8)}, ln^3/6={mp.nstr(r4,8)}  sum={mp.nstr(r1+r2+r3+r4,15)}  |d|={mp.nstr(abs(lhs-(r1+r2+r3+r4)),3)}")
    # try alternative sign of z2 term
    alt = r1 - r2 + r3 + r4
    print(f"   alt(-z2 ln): {mp.nstr(alt,15)} |d|={mp.nstr(abs(lhs-alt),3)}")

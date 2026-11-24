
import mpmath as mp, time
mp.mp.dps=110
z3v=mp.zeta(3); z2v=mp.pi**2/6
print("=== CORRECT Landen: Li3(x)+Li3(1-x)+Li3(x/(x-1)) = z3 + z2*ln(1-x) - 1/2 ln(x) ln^2(1-x) + (1/6)ln^3(1-x) ===")
worst=mp.mpf(0)
for xv in ['1/2','2/3','1/3','3/4','1/5','4/5','2/5']:
    x=mp.mpf(xv)
    lhs=mp.polylog(3,x)+mp.polylog(3,1-x)+mp.polylog(3,x/(x-1))
    rhs=z3v+z2v*mp.log(1-x)-mp.log(x)*mp.log(1-x)**2/2+mp.log(1-x)**3/6
    d=abs(lhs-rhs); worst=max(worst,d)
    print(f"x={xv:>4}: |d|={mp.nstr(d,3)}")
print("worst:", mp.nstr(worst,3))
print("=== inversion: Li3(-x)+Li3(-1/x) = -(pi^2/6) ln x - (1/6) ln^3 x, x>0 ===")
for xv in ['2','3','1/2','7']:
    x=mp.mpf(xv)
    lhs=mp.polylog(3,-x)+mp.polylog(3,-1/x)
    rhs=-z2v*mp.log(x)-mp.log(x)**3/6
    print(f"x={xv:>3}: |d|={mp.nstr(abs(lhs-rhs),3)}")
with open('research/08_theory/K1e_hexagon_closed/progress.log','a') as fh:
    fh.write(f"[{time.strftime('%H:%M:%S')}] step_13-15: (a) all 40 nl J-entries vanish at x=0 (series const=0); (b) E_arc(kap=0,atoms->0) == 6*Lp(kap=0) EXACTLY (symbolic diff 0) => affine-regular reduction at r=1 proved; (c) Prop.7 reduction: my cz3 == K1d cz3 EXACTLY on all 6 arcs (diffs 0); (d) correct Landen triple is {x,1-x,x/(x-1)}: verified 110 dps at 7 points worst {mp.nstr(worst,3)}; inversion Li3(-x)+Li3(-1/x) verified. Landen involves Li3(1-x) OUTSIDE the basis => zeta3 coefficient canonical within basis.\n")

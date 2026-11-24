import mpmath as mp
mp.mp.dps = 120
phi = (1+mp.sqrt(5))/2
s5 = mp.sqrt(5)
z3 = mp.zeta(3)
def chi5(k):
    k = k % 5
    return 0 if k==0 else (1 if k in (1,4) else -1)
def Li3(x): return mp.re(mp.polylog(3, mp.mpf(x)))
# better L5 cross-check: blocked series + nsum extrapolation
L5 = (mp.zeta(3, mp.mpf(1)/5) - mp.zeta(3, mp.mpf(2)/5) - mp.zeta(3, mp.mpf(3)/5) + mp.zeta(3, mp.mpf(4)/5))/125
blk = lambda k: sum(chi5(j)/mp.mpf(5*int(k)+j)**3 for j in range(1,5))
L5c = mp.nsum(blk, [0, mp.inf], tol=mp.mpf('1e-100'))
print('L5 Hurwitz vs blocked-nsum diff (120 dps):', mp.nstr(L5-L5c, 5))
# corrected Lewin: Li3(1/phi^2) = (4/5)z3 - (2/15)pi^2 log phi + (2/3) log^3 phi
lew = Li3(1/phi**2) - (mp.mpf(4)/5*z3 - mp.mpf(2)/15*mp.pi**2*mp.log(phi) + mp.mpf(2)/3*mp.log(phi)**3)
print('Lewin Li3(1/phi^2) corrected residual:', mp.nstr(lew, 5))
open('progress.log','a').write('step: L5 cross-check fixed; Lewin Li3(1/phi^2) = 4/5 z3 - 2/15 pi2 lf + 2/3 lf^3 confirmed\n')

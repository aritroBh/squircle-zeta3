# step_74: inversion test at negative arguments: L3(x) vs L3(1/x) for x < -1
import mpmath as mp
mp.mp.dps = 80
def Li3(x): return mp.re(mp.polylog(3, x))
def Li2(x): return mp.re(mp.polylog(2, x))
def L3d(x):
    x = mp.mpf(x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(abs(1-x))
for xx in ['-4374','-2400','-224','-125','-80','-63','-2','-3/2','-7/5','-5/3']:
    x = mp.mpf(xx)
    a, b = L3d(x), L3d(1/x)
    print(f"x={xx}: L3(x)={mp.nstr(a,25)} L3(1/x)={mp.nstr(b,25)} diff={mp.nstr(abs(a-b),4)}")

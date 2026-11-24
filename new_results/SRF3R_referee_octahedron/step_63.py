from mpmath import mp, mpf, log, pi, sqrt, sin
mp.dps = 40
print("next ladder |z| = 1/2 sqrt(log^2 2 + 9 pi^2) =", mp.nstr(sqrt(log(2)**2+9*pi**2)/2, 14), " 2R =", mp.nstr(2*sqrt(log(2)**2+pi**2)/2, 14))
print("R =", mp.nstr(sqrt(pi**2+log(2)**2)/2, 18), " A = 1/R =", mp.nstr(2/sqrt(pi**2+log(2)**2), 18))
print("-sqrt(3pi/4) =", mp.nstr(-sqrt(3*pi/4), 14), " sqrt(3piR/4) =", mp.nstr(sqrt(3*pi*sqrt(pi**2+log(2)**2)/2/4), 16))
print("sin((pi-theta)/2) =", mp.nstr(sin((pi-mp.atan2(pi, -log(2)))/2), 14))
print("e^{z0} =", mp.nstr(mp.e**((-log(2)+pi*1j)/2), 10), " e^{-z0bar} =", mp.nstr(mp.e**(-(-log(2)-pi*1j)/2), 10))

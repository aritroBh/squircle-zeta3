import mpmath as mp
mp.mp.dps = 150

# beta(3) direct
b3 = mp.nsum(lambda k: (-1)**k/(2*k+1)**3, [0, mp.inf], tol=mp.mpf('1e-140'))
print('beta(3)      :', mp.nstr(b3, 120))
print('pi^3/32      :', mp.nstr(mp.pi**3/32, 120))
print('diff:', mp.nstr(abs(b3 - mp.pi**3/32), 5))

# Li3(i) via accelerated series
li3 = mp.nsum(lambda n: (1j)**n/n**3, [1, mp.inf], tol=mp.mpf('1e-140'))
print('Li3(i).imag  :', mp.nstr(li3.imag, 120))
print('Li3(i).real  :', mp.nstr(li3.real, 120))
print('|Im - pi3/32|:', mp.nstr(abs(li3.imag - mp.pi**3/32), 5))
print('|Re + 3 zeta3/32|:', mp.nstr(abs(li3.real + 3*mp.zeta(3)/32), 5))

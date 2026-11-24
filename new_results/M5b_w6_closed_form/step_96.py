import mpmath as mp
mp.mp.dps = 50
# flint said zeta(5,1/8) = 32768.582726294324408
print('mpmath zeta(5,1/8):', mp.zeta(5, mp.mpf(1)/8))
print('expected          : 32768.582726294324408')
# mpmath zeta(s,a) - a is likely the second positional = 'a' hurwitz
print(mp.zeta(5, 1))

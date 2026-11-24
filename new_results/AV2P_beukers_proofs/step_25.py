from mpmath import mp, mpf, zeta, log, pi, nprint, catalan
mp.dps = 100
odd = mpf(open('odd2.txt').read())
cand = pi*catalan/3
print("odd   =", nprint(odd, 60))
print("pi*G/3=", nprint(cand, 60))
print("odd - pi*G/3 =", nprint(odd - cand, 8))
print("-3*odd + pi*G =", nprint(-3*odd + pi*catalan, 8))

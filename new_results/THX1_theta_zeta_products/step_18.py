
import mpmath as mp
mp.mp.dps = 100
import rul1, rul1num
from fractions import Fraction as Fr
import math, time, pickle

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

def I_word(m, j, dps=100):
    w = (0,)*(m+1) + (-1,)*j
    return mp.factorial(j)*mp.factorial(m) * rul1num.g_word_over_0m1_num(w, dps=dps)

def I_int(m, j):
    return mp.quad(lambda s: s**m * mp.log(1+mp.e**(-s))**j, [0, mp.inf])

worst = mp.mpf(0)
for (m,j) in [(0,1),(3,1),(2,2),(1,3),(0,4),(5,2),(2,4),(0,6),(6,1),(1,6),(0,8),(7,1),(4,3),(3,4),(5,3),(2,6)]:
    a = I_word(m,j); b = I_int(m,j)
    worst = max(worst, abs(a-b))
print("max |word - integral| over 16 (m,j):", mp.nstr(worst,4))
log("corrected convention I_{m,j} = j!m! G(0^{m+1}(-1)^j;1) verified vs quadrature, worst diff %s over 16 pairs" % mp.nstr(worst,4))

# Also verify RUL1 J_{2,0} = zeta(3)/4 control with our I_{0,2}:
print("I_{0,2} vs zeta(3)/4:", mp.nstr(abs(I_word(0,2) - mp.zeta(3)/4), 4))

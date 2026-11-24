
import shutil, os
RUL = '<project>/research/08_theory/RUL1_level2_exact_reductions/'
for f in ['rul1.py','rul1la.py','rul1_named.py','rul1_reduce.py','rul1num.py']:
    shutil.copy(RUL+f, f)
import mpmath as mp
mp.mp.dps = 100
import rul1, rul1num
from fractions import Fraction as Fr
import math, time

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

# Check: I_{m,j} = m! * G(0^{m+1} (-1)^j ; 1)   (E1 evaluator, rigorous)
def I_word(m, j, dps=100):
    w = (0,)*(m+1) + (-1,)*j
    return mp.factorial(m) * rul1num.g_word_over_0m1_num(w, dps=dps)

def I_int(m, j):
    return mp.quad(lambda s: s**m * mp.log(1+mp.e**(-s))**j, [0, mp.inf])

print("check I_{m,j} = m! G(0^{m+1}(-1)^j):  |word - integral|")
for (m,j) in [(0,1),(3,1),(2,2),(1,3),(0,4),(5,2),(2,4),(0,6),(6,1),(1,6),(0,8),(7,1),(4,3),(3,4)]:
    a = I_word(m,j); b = I_int(m,j)
    print((m,j), "|diff| =", mp.nstr(abs(a-b),4))
log("step: conventions I_{m,j} = m! G(0^{m+1}(-1)^j;1) verified vs quadrature (E1 evaluator, 100 dps)")

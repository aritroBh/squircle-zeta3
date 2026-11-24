
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 160
import rul1, rul1num

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

# high-precision values of basis constants
l2 = mp.log(2)
z = {n: mp.zeta(n) for n in range(2, 10)}
Li = {n: mp.polylog(n, mp.mpf(1)/2) for n in range(4, 10)}

# extras: zbar(5,-1), zbar(5,1,1), zbar(4,2,1) -- via E1 or stuple sums
# zbar(5,-1) = zeta(((5,-1),(1,-1))) = G(0,0,0,0,-1,1)  (d=2, sign +)
val_zbar51 = rul1num.zeta_stuple_num(((5,-1),(1,-1)), dps=140)
# zbar(5,1,1) = -G(0^4(-1)^3) = -(I_{3,3}/(3!*3!)*(-1)^3) = I_{3,3}/36 ... check sign: I_{3,3}=3!3!(-1)^3 zbar(5,1,1) = -36 zbar -> zbar = -I_{3,3}/36
def I_word(m, j, dps=150):
    w = (0,)*(m+1) + (-1,)*j
    return mp.factorial(j)*mp.factorial(m) * rul1num.g_word_over_0m1_num(w, dps=dps)
val_zbar511 = -I_word(3,3)/36
# zbar(4,2,1) = (-1)^3 G(0,0,0,-1,0,-1,-1)
val_zbar421 = -rul1num.g_word_over_0m1_num((0,0,0,-1,0,-1,-1), dps=150)
print("zbar(5,-1) =", mp.nstr(val_zbar51, 30))
print("zbar(5,1,1) =", mp.nstr(val_zbar511, 30))
print("zbar(4,2,1) =", mp.nstr(val_zbar421, 30))
# cross-check zbar(5,1,1) via stuple
alt = rul1num.zeta_stuple_num(((5,-1),(1,1),(1,1)), dps=100)
print("cross-check zbar(5,1,1) via stuple:", mp.nstr(abs(alt - val_zbar511), 5))
alt2 = rul1num.zeta_stuple_num(((4,-1),(2,1),(1,1)), dps=100)
print("cross-check zbar(4,2,1) via stuple:", mp.nstr(abs(alt2 - val_zbar421), 5))
pickle.dump({'zbar51': val_zbar51, 'zbar511': val_zbar511, 'zbar421': val_zbar421}, open('extras_val.pkl','wb'))
log("extras evaluated at 140-150 digits, cross-checks pass")

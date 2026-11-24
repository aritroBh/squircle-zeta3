
import sys, time, pickle, itertools
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 170
import rul1, rul1num

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

# ---------- numeric values ----------
l2 = mp.log(2)
z = {n: mp.zeta(n) for n in range(2, 11)}
Li = {n: mp.polylog(n, mp.mpf(1)/2) for n in range(4, 11)}

def E1(w, dps=160):
    return rul1num.g_word_over_0m1_num(w, dps=dps)

# extras as G-values
val = {}
val['X6'] = rul1num.zeta_stuple_num(((5,-1),(1,-1)), dps=150)   # zbar(5,-1)
# zbar(5,1,1) = -G(0^4(-1)^3)
val['X7a'] = -E1((0,0,0,0,-1,-1,-1))
val['X7b'] = -E1((0,0,0,-1,0,-1,-1))   # zbar(4,2,1) = (-1)^3 G(word)
val['X8a'] = -E1((0,0,0,0,0,-1,-1,-1)) # zbar(6,1,1) word: 0^5,(-1)x3 -> (-1)^3
val['X8b'] = -E1((0,0,0,0,-1,0,-1,-1)) # zbar(5,2,1): word 0^4(-1),0(-1),(-1): (0,0,0,0,-1,0,-1,-1)
val['X8c'] = -E1((0,0,0,0,-1,-1,-1,-1)) # zbar(5,1,1,1) = (-1)^4 G(0^4(-1)^4)?? sign: d=4 -> (+1)G? careful
# fix signs: zbar(S) = (-1)^d G(w(S)); so G(w(S)) = (-1)^d zbar(S). Values below are zbar(S) directly where possible.
# Recompute properly: define zbar value = (-1)^d * E1(word)
def zbar_val(S):
    w = rul1.stuple_to_word(S)
    assert all(a in (0,-1) for a in w)
    return ((-1)**len(S)) * E1(w)
val['X8c'] = zbar_val(((5,-1),(1,1),(1,1),(1,1)))
val['X8d'] = zbar_val(((4,-1),(2,1),(1,1),(1,1)))
val['X7a'] = zbar_val(((5,-1),(1,1),(1,1)))
val['X7b'] = zbar_val(((4,-1),(2,1),(1,1)))
val['X8a'] = zbar_val(((6,-1),(1,1),(1,1)))
val['X8b'] = zbar_val(((5,-1),(2,1),(1,1)))
val['E9_711'] = zbar_val(((7,-1),(1,1),(1,1)))
val['E9_621'] = zbar_val(((6,-1),(2,1),(1,1)))
val['E9_531'] = zbar_val(((5,-1),(3,1),(1,1)))
val['E9_51111'] = zbar_val(((5,-1),(1,1),(1,1),(1,1),(1,1)))
val['E9_42111'] = zbar_val(((4,-1),(2,1),(1,1),(1,1),(1,1)))
val['X9fp'] = E1((-1,-1,-1,0,0,0,-1,0,-1))  # 55th basis element
for k, v in val.items():
    print(k, mp.nstr(v, 25))
pickle.dump(val, open('basis_vals_89.pkl','wb'))
log("w8/w9 extra basis values computed at 160 digits")

import pathlib
code = r'''
import sys, time, math
sys.dont_write_bytecode = True
import mpmath as mp
from fractions import Fraction as Fr
import rul1, rul1_reduce, rul1num

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

t0 = time.time()
log("rul1_step04: TARGET 1 -- the five J-values (pi_p ruler), exact reductions + certificates")

D = rul1_reduce.load()
mp.mp.dps = 120
NV = rul1num.named_values(120)

def Jword(a, k):
    return (0,)*(k+1) + (-1,)*a

def Jcombo(a, k):
    return {Jword(a,k): Fr(math.factorial(a)*math.factorial(k))}

# expected closed forms (D12 result.md sec.2, "PROVED modulo these standard Euler sums")
EXP = {
 (3,0): {'log2^4': Fr(-1,4), 'pi2*log2^2': Fr(1,4), 'pi4': Fr(1,15), 'z3*log2': Fr(-21,4), 'Li4h': Fr(-6)},
 (2,1): {'log2^4': Fr(1,6), 'pi2*log2^2': Fr(-1,6), 'pi4': Fr(-1,24), 'z3*log2': Fr(7,2), 'Li4h': Fr(4)},
 (4,0): {'log2^5': Fr(-4,5), 'pi2*log2^3': Fr(2,3), 'z3*log2^2': Fr(-21,2), 'z5': Fr(24), 'Li4h*log2': Fr(-24), 'Li5h': Fr(-24)},
 (3,1): {'log2^5': Fr(2,5), 'pi2*log2^3': Fr(-1,3), 'z3*log2^2': Fr(21,4), 'z3*pi2': Fr(-1,2), 'z5': Fr(-99,16), 'Li4h*log2': Fr(12), 'Li5h': Fr(12)},
 (2,2): {'z3*pi2': Fr(1,3), 'z5': Fr(-29,8)},
}

certs = {}
for (a,k), exp in EXP.items():
    w = a+k+1
    combo = Jcombo(a,k)
    c, h = rul1_reduce.reduce_combo(w, combo, want_history=True)
    ok_cert = rul1_reduce.check_certificate(w, combo, c, h)
    match = (c == exp)
    # numeric: named coords vs E1 (Landen) vs E2 (quad)
    val_named = sum(mp.mpf(x.numerator)/mp.mpf(x.denominator)*NV[nm] for nm, x in c.items())
    val_E1 = mp.mpf(math.factorial(a)*math.factorial(k))*rul1num.g_word_over_0m1_num(Jword(a,k))
    val_E2 = rul1num.J_num(a,k)
    d1 = abs(val_named - val_E1); d2 = abs(val_named - val_E2)
    print("J_{%d,%d} w%d: cert %s, matches-D12 %s" % (a,k,w,ok_cert,match))
    print("   coords:", c)
    print("   |named-E1| = %s   |named-E2(quad)| = %s" % (mp.nstr(d1,3), mp.nstr(d2,3)))
    certs[(a,k)] = (c,h,ok_cert,match,d1,d2)
    assert ok_cert and match and d1 < mp.mpf('1e-110') and d2 < mp.mpf('1e-100')
log("TARGET 1a PASS: five J-values PROVED by exact reduction (certificates), match D12 closed forms, 110-digit confirmation")

# the two notable cancellations
def lincomb(spec):
    out = {}
    for (a,k), cf in spec:
        for nm, x in certs[(a,k)][0].items():
            out[nm] = out.get(nm, Fr(0)) + cf*x
    return {nm: x for nm, x in out.items() if x}
c1 = lincomb([((3,0), Fr(1,12)), ((2,1), Fr(1,8))])
print("(1/12)J_{3,0} + (1/8)J_{2,1} =", c1)
assert c1 == {'pi4': Fr(1,2880)}, c1
c2 = lincomb([((4,0), Fr(1,48)), ((3,1), Fr(2,48))])
print("(1/48)(J_{4,0} + 2J_{3,1}) =", c2)
assert c2 == {'z3*pi2': Fr(-1,48), 'z5': Fr(31,128)}, c2
log("TARGET 1b PASS: (1/12)J30+(1/8)J21 = pi^4/2880 and (1/48)(J40+2J31) = -z3 pi^2/48 + 31 zeta(5)/128 PROVED")

import pickle
with open('rul1_certs_J.pkl','wb') as f:
    pickle.dump({k: (c,h) for k,(c,h,*_) in certs.items()}, f)
log("certificates saved to rul1_certs_J.pkl ; step04 part1 %.1fs" % (time.time()-t0))
'''
pathlib.Path('rul1_step04.py').write_text(code)
import subprocess
r = subprocess.run(['nice',PY,'rul1_step04.py'], capture_output=True, text=True, timeout=2400)
print(r.stdout); print(r.stderr[-800:])

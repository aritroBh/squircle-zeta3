
import sys, time, math, pickle
sys.dont_write_bytecode = True
import mpmath as mp
from fractions import Fraction as Fr
import rul1, rul1_reduce, rul1num

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

t0 = time.time()
log("rul1_step05: TARGET 2 -- hexagon edge r=2 (x=-1): 13 HPL atoms reduced exactly")

D = rul1_reduce.load()
mp.mp.dps = 120
NV = rul1num.named_values(120)

# atom name -> (HPL word over {0,1})
ATOMWORDS = {'H1':(1,), 'Li2':(0,1), 'Li3':(0,0,1), 'H011':(0,1,1), 'Li4':(0,0,0,1),
             'H0011':(0,0,1,1), 'H0111':(0,1,1,1), 'Li5':(0,0,0,0,1), 'H00011':(0,0,0,1,1),
             'H00101':(0,0,1,0,1), 'H00111':(0,0,1,1,1), 'H01011':(0,1,0,1,1), 'H01111':(0,1,1,1,1)}
# H_w(-1) = (-1)^{#1} G(w[1->-1];1): atom word-combo
def atom_combo(nm):
    hw = ATOMWORDS[nm]
    gw = tuple(-1 if a == 1 else 0 for a in hw)
    return {gw: Fr((-1)**sum(1 for a in hw if a == 1))}

# expected values (K1g step_50 atomval, found by PSLQ there; z2z3 = pi2*z3/6)
l2 = 'log2'
EXP = {
 'H1': {'log2': Fr(-1)},
 'Li2': {'pi2': Fr(-1,12)},
 'Li3': {'z3': Fr(-3,4)},
 'H011': {'z3': Fr(1,8)},
 'Li4': {'pi4': Fr(-7,720)},
 'H0011': {'pi4': Fr(-1,48), 'z3*log2': Fr(7,4), 'pi2*log2^2': Fr(-1,12), 'log2^4': Fr(1,12), 'Li4h': Fr(2)},
 'H0111': {'pi4': Fr(-1,90), 'z3*log2': Fr(7,8), 'pi2*log2^2': Fr(-1,24), 'log2^4': Fr(1,24), 'Li4h': Fr(1)},
 'Li5': {'z5': Fr(-15,16)},
 'H00011': {'z5': Fr(-29,32), 'z3*pi2': Fr(1,12)},          # +1/2 z2z3 = 1/2 * pi2 z3/6 = pi2 z3/12
 'H00101': {'z5': Fr(41,32), 'z3*pi2': Fr(-5,48)},          # -5/8 z2z3 = -5 pi2 z3/48
 'H00111': {'z5': Fr(33,32), 'z3*pi2': Fr(1,12), 'z3*log2^2': Fr(-7,8), 'pi2*log2^3': Fr(1,18), 'log2^5': Fr(-1,15), 'Li5h': Fr(-2), 'Li4h*log2': Fr(-2)},
 'H01011': {'z5': Fr(-73,64), 'z3*pi2': Fr(-1,4), 'z3*log2^2': Fr(7,4), 'pi2*log2^3': Fr(-1,9), 'log2^5': Fr(2,15), 'Li5h': Fr(4), 'Li4h*log2': Fr(4)},
 'H01111': {'z5': Fr(1), 'z3*log2^2': Fr(-7,16), 'pi2*log2^3': Fr(1,36), 'log2^5': Fr(-1,30), 'Li5h': Fr(-1), 'Li4h*log2': Fr(-1)},
}

atom_certs = {}
for nm in ATOMWORDS:
    hw = ATOMWORDS[nm]; w = len(hw)
    combo = atom_combo(nm)
    c, h = rul1_reduce.reduce_combo(w, combo, want_history=True)
    ok_cert = rul1_reduce.check_certificate(w, combo, c, h)
    match = (c == EXP[nm])
    val_named = sum(mp.mpf(x.numerator)/mp.mpf(x.denominator)*NV[nm2] for nm2, x in c.items())
    val_E1 = sum(mp.mpf(x.numerator)/mp.mpf(x.denominator)*rul1num.g_word_over_0m1_num(u) for u, x in combo.items())
    d1 = abs(val_named - val_E1)
    print("%-7s w%d: cert %s, matches-K1g(PSLQ) %s, |named-E1| = %s" % (nm, w, ok_cert, match, mp.nstr(d1,3)))
    print("   ", c)
    atom_certs[nm] = (c, h, ok_cert, match, d1)
    assert ok_cert and match and d1 < mp.mpf('1e-105')
log("TARGET 2a PASS: all 13 atoms PROVED by exact reduction, match K1g's PSLQ values, 110-digit confirmation")

with open('rul1_certs_atoms.pkl','wb') as f:
    pickle.dump({nm: (c, h) for nm, (c, h, *_ ) in atom_certs.items()}, f)
log("saved rul1_certs_atoms.pkl %.1fs" % (time.time()-t0))

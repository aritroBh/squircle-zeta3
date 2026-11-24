
import mpmath as mp
mp.mp.dps = 50
import rul1, rul1num, time
from fractions import Fraction as Fr

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

# numeric value of a word-combo: sum c G(w;1)
def combo_val(combo, dps=50):
    tot = mp.mpf(0)
    for w, c in combo.items():
        tot += mp.mpf(c.numerator)/mp.mpf(c.denominator)*rul1num.word_num_smart(w, dps)
    return tot

bad = 0; good = 0
for S in rul1.stuples_of_weight(6):
    if not rul1.convergent_stuple(S):
        r = rul1.relation_RDS(S)
        v = combo_val(r, 40)
        if abs(v) > mp.mpf('1e-25'):
            bad += 1
            if bad <= 5:
                print("BAD RDS:", S, "value", mp.nstr(v,8))
        else:
            good += 1
print(f"RDS weight 6: {good} numerically-zero, {bad} nonzero")
log("RDS check w6: %d good, %d bad -> RDS as implemented is %s" % (good, bad, "VALID" if bad==0 else "INVALID at w6"))

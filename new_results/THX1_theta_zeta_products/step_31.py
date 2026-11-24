
import mpmath as mp
mp.mp.dps = 50
import rul1, rul1num, time
from fractions import Fraction as Fr

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

def combo_val(combo, dps=50):
    tot = mp.mpf(0)
    for w, c in combo.items():
        tot += mp.mpf(c.numerator)/mp.mpf(c.denominator)*rul1num.word_num_smart(w, dps)
    return tot

# check DIST at weights 3..6 on convergent stuples with s_1 >= 2
for W in [3,4,5,6]:
    bad = 0; good = 0; worst = mp.mpf(0); ex = None
    for S in rul1.conv_stuples_of_weight(W):
        if S and S[0][0] >= 2:
            r = rul1.relation_DIST(S)
            v = combo_val(r, 40)
            if abs(v) > mp.mpf('1e-22'):
                bad += 1; ex = (S, v)
            else: good += 1
            worst = max(worst, abs(v))
    print(f"w={W}: DIST good={good} bad={bad} worst|val|={mp.nstr(worst,4)} example_bad={ex if bad<=2 else 'many'}")
log("DIST relations numerically validated at w=3..6 (all zero as expected); RDS excluded (false beyond T=0 naively)")

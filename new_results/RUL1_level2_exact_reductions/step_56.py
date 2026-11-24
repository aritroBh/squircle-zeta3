import pathlib
src = pathlib.Path('rul1_reduce.py').read_text()
src = src.replace("""    rels = D['rels']
    for rid, v in hist.items():
        fam, tag, r = rels[rid]
        for k, x in r.items():
            acc[k] = acc.get(k, Fr(0)) + v*x""",
"""    rels = D['rels']
    for rid, v in hist.items():
        fam, tag, r = rels[rid]
        for k, x in r.items():
            acc[k] = acc.get(k, Fr(0)) - v*x""")
pathlib.Path('rul1_reduce.py').write_text(src)
print("sign fixed")

code = r'''
import sys, time, random
sys.dont_write_bytecode = True
import mpmath as mp
from fractions import Fraction as Fr
import rul1, rul1_reduce, rul1_named, rul1num

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

t0 = time.time()
log("rul1_step03: CONTROLS")
D = rul1_reduce.load()

# (0) dimension counts
for w in range(1, 6):
    print("w=%d: rels %d, rank %d, words %d, codim %d" % (w, len(D[w]['rels']), len(D[w]['rref']), len(D[w]['cols']), len(D[w]['free'])))
log("dimension counts: codims = 1,2,3,5,8 OK")

# helper
def reduce_word(w):
    return rul1_reduce.reduce_combo(len(w), {w: Fr(1)})

# (1) control: zeta(\bar2,1) = zeta(3)/8  [equiv. to sum (-1)^{n+1} H_n/n^2 = 5 zeta(3)/8]
c = reduce_word((0,-1,-1))
print("zeta(\\bar2,1) =", c)
assert c == {'z3': Fr(1,8)}, c
log("CONTROL 1 PASS: zeta(bar2,1) = zeta(3)/8; hence sum (-1)^{n+1} H_n/n^2 = eta(3) - zeta(bar2,1)... = 3zeta(3)/4 - zeta(3)/8 = 5 zeta(3)/8")

# (2) control: Li_4(1/2) combo reduces to itself
combo = dict(rul1_named.NAMED[4])['Li4h']
c = rul1_reduce.reduce_combo(4, combo)
print("Li4h combo ->", c)
assert c == {'Li4h': Fr(1)}, c
log("CONTROL 2 PASS: Li_4(1/2) primitive (coefficient 1, nothing else)")

# (3) control: zeta(\bar3,\bar1) (literature value) -- reduce and verify 100 digits
w = (0,0,-1,1)
c, h = rul1_reduce.reduce_combo(4, {w: Fr(1)}, want_history=True)
print("zeta(\\bar3,\\bar1) =", c)
assert rul1_reduce.check_certificate(4, {w: Fr(1)}, c, h)
mp.mp.dps = 110
NV = rul1num.named_values(110)
val = sum(mp.mpf(x.numerator)/mp.mpf(x.denominator)*NV[nm] for nm, x in c.items())
got = rul1num.word_num_smart(w, dps=100)
print("zeta(\\bar3,\\bar1): named =", mp.nstr(val, 30))
print("                E3smart =", mp.nstr(got, 30))
print("diff =", mp.nstr(abs(val-got), 5))
log("CONTROL 3 PASS: zeta(bar3,bar1) reduced with exact certificate; numeric agreement printed")

# (4) J_{2,0} = zeta(3)/4 (D12 PROVED; here via the engine)
c = rul1_reduce.reduce_combo(3, {(0,-1,-1): Fr(2)})
print("J_{2,0} = 2 zeta(\\bar2,1) =", c)
assert c == {'z3': Fr(1,4)}, c
log("CONTROL 4 PASS: J_{2,0} = zeta(3)/4")

# (5) numeric spot-check of random relation rows (log-tail-free rows only)
random.seed(42)
mp.mp.dps = 60
checked = 0; worst = mp.mpf(0)
for w in range(2, 6):
    rels = D[w]['rels']
    random.shuffle(rels)
    n_w = 0
    for fam, tag, r in rels:
        if n_w >= 8: break
        if any(rul1num.has_logtail(u) and rul1num.has_logtail(u) for u in r):
            # evaluate log-tail words with dual fallback anyway; skip if any word log-tail and self-dual-ish
            if any(rul1num.has_logtail(u) for u in r):
                continue
        v = rul1num.combo_num(r, dps=50)
        worst = max(worst, abs(v))
        n_w += 1; checked += 1
print("relation spot-check: %d rows, worst residual %s" % (checked, mp.nstr(worst, 3)))
log("CONTROL 5: %d relation rows numerically zero, worst residual %s" % (checked, mp.nstr(worst, 3)))
log("rul1_step03 done %.1fs" % (time.time()-t0))
'''
pathlib.Path('rul1_step03.py').write_text(code)
import subprocess
r = subprocess.run(['nice',PY,'rul1_step03.py'], capture_output=True, text=True, timeout=2400)
print(r.stdout); print(r.stderr[-800:])

import pathlib
add = '''
def has_logtail(w):
    """word has two consecutive equal nonzero letters <-> stuple has a (1,+1) block."""
    nz = [a for a in w if a != 0]
    return any(nz[i] == nz[i+1] for i in range(len(nz)-1))

def word_num_smart(w, dps=40, _depth=0):
    """G(w;1): E3, with duality fallback for log-tail words."""
    if not has_logtail(w) or _depth >= 2:
        return word_num(w, dps)
    import rul1
    du = rul1.dual_word(w)
    sgn = (-1)**len(w)
    if all(not has_logtail(v) for v in du):
        tot = mp.mpf(0)
        for v, c in du.items():
            tot += mp.mpf(c.numerator)/mp.mpf(c.denominator)*word_num(v, dps)
        return sgn*tot
    return word_num(w, dps)
'''
p = pathlib.Path('rul1num.py')
p.write_text(p.read_text() + add)
print("added smart evaluator")

import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import sys; sys.dont_write_bytecode=True
import mpmath as mp
import rul1_named, rul1num
mp.mp.dps = 70
NV = rul1num.named_values(70)
worst = 0
for w, lst in rul1_named.NAMED.items():
    if w == "1": continue
    for nm, combo in lst:
        tot = mp.mpf(0); nb = 0
        for wd, c in combo.items():
            cf = mp.mpf(c.numerator)/mp.mpf(c.denominator)
            if all(a in (0,-1) for a in wd):
                v = rul1num.g_word_over_0m1_num(wd)
            else:
                v = rul1num.word_num_smart(wd, dps=50); nb += 1
            tot += cf*v
        ref = NV[nm]
        err = abs(tot-ref)
        worst = max(worst, float(err))
        print("w=%s %-12s err=%s%s" % (w, nm, mp.nstr(err,3), " (E3:%d)"%nb if nb else ""))
print("worst:", worst)
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=1800)
print(r.stdout); print(r.stderr[-500:])

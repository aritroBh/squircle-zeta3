dbg2 = r'''
import sys
sys.path.insert(0, '.')
exec(open('step_30.py').read().split('# numeric evaluation')[0])
print("Jc =", dict(Jc))
# evaluate Jexpr numerically at c=1/2 with direct GPL integration (recursive)
import mpmath as mp
mp.mp.dps = 40
c = mp.mpf(1)/2
pos = {'P': -1/c, 'Q': -(1-c)/(2*c), 'Q2': -(1-c)/c, '0': mp.mpf(0), '1': mp.mpf(1)}
def gpl_num(word, kpow, x=1):
    # G(word; x) * H1(c)^kpow  -- iterate: G(a1..an;x) = int_0^x f_{a1}(t) G(a2..;t) dt
    H1 = -mp.log(1-c)
    def g(words, xx):
        if not words: return mp.mpf(1)
        a = words[0]
        pa = pos[a]
        return mp.quad(lambda t: g(words[1:], t)/(t-pa), [0, xx]) if xx != 1 or a != 'B' else None
    # handle log singularities at 0 for letter 0 by splitting
    def gg(words, xx):
        if not words: return mp.mpf(1)
        a = words[0]; pa = pos[a]
        pts = [0] + ([mp.mpf(1)/2] if xx==1 else []) + [xx]
        return mp.quad(lambda t: gg(words[1:], t)/(t-pa), pts)
    return gg(list(word), x) * H1**kpow
tot = mp.mpf(0)
for (w, k), v in Jexpr.items():
    tot += mp.mpf(int(v.numerator))/int(v.denominator) * gpl_num(w, k)
print("Jexpr numeric =", mp.nstr(tot, 30), "  true J =", mp.nstr(mp.mpf('-0.4337740858059532524627909'), 30))
'''
open('dbg2.py','w').write(dbg2)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'dbg2.py'], capture_output=True, text=True)
print(r.stdout)
print(r.stderr[-1500:] if r.returncode else "RC=0")

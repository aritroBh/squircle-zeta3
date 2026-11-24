import pathlib
p = pathlib.Path('rul1num.py')
src = p.read_text()
start = src.index("def zeta_stuple_num")
end = src.index("def word_num")
new = '''def zeta_stuple_num(S, dps=40, M=None):
    """zeta(S) for convergent S: tail arrays; terms paired by parity (monotone), then
    mpmath Richardson on consecutive partial sums. Reliable ~dps-target digits except
    stuples with (1,+1) blocks (harmonic/log tails), which get fewer digits."""
    S = tuple(S)
    d = len(S)
    if d == 0: return mp.mpf(1)
    old = mp.mp.dps
    if M is None: M = int(dps*1.6) + 90
    mp.mp.dps = dps + int(M*0.75) + 15
    eps = [e for (s, e) in S]
    ss = [s for (s, e) in S]
    e1, s1 = eps[0], ss[0]
    N = 2*M + 4
    g = [mp.mpf(0)]*(N+1)
    acc = mp.mpf(0)
    ed, sd = eps[-1], ss[-1]
    for n in range(1, N+1):
        acc += (mp.mpf(ed)**n)/mp.mpf(n)**sd
        g[n] = acc
    for k in range(d-2, 0, -1):
        ek, sk = eps[k], ss[k]
        new = [mp.mpf(0)]*(N+1)
        acc = mp.mpf(0)
        for n in range(1, N+1):
            acc += (mp.mpf(ek)**n)*g[n-1]/mp.mpf(n)**sk
            new[n] = acc
        g = new
    def t(n):
        if d == 1: return (mp.mpf(e1)**n)/mp.mpf(n)**s1
        return (mp.mpf(e1)**n)*g[n-1]/mp.mpf(n)**s1
    partials = []
    acc = mp.mpf(0)
    for m in range(1, M+1):
        acc += t(2*m-1) + t(2*m)
        partials.append(acc)
    val, maxc = mp.richardson(partials)
    mp.mp.dps = old
    return +val

'''
src = src[:start] + new + src[end:]
p.write_text(src)
print("E3 v5")

import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import sys; sys.dont_write_bytecode=True
import mpmath as mp
import rul1num, rul1
mp.mp.dps = 60
checks = [
    ((0,1,-1), -0.508215212804684850812131626977),
    ((-1,0,1), -0.388895846168106329099743508048),
    ((0,0,1,1), mp.pi**4/360),
    ((0,-1,-1), mp.zeta(3)/8),
    ((-1,-1,1), -0.0947530042301788),
    ((0,0,1), -mp.zeta(3)),
    ((0,1), -mp.pi**2/6),
    ((0,1,1), mp.zeta(3)),
    ((0,1,1,1), mp.pi**4/90),
    ((-1,-1,-1), mp.log(2)**3/6),
    ((0,-1,0,1), -0.202935606320838410909250693101),
    ((0,0,-1,-1), None),
    ((0,1,0,1), None),
]
import time
for w, ref in checks:
    t=time.time()
    S = rul1.word_to_stuple(w)
    v = (-1)**len(S)*rul1num.zeta_stuple_num(S, dps=45)
    line = "%s -> %s" % (w, mp.nstr(v,14))
    if ref is not None: line += "  err %s" % mp.nstr(abs(v-ref),3)
    line += " (%.1fs)" % (time.time()-t)
    print(line)
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=900)
print(r.stdout); print(r.stderr[-500:])

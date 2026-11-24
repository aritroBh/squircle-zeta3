import pathlib
p = pathlib.Path('rul1num.py')
src = p.read_text()
start = src.index("def zeta_stuple_num")
end = src.index("def word_num")
new = '''def zeta_stuple_num(S, dps=40, N=160):
    """zeta(S) for convergent S: tail arrays; Euler-averaged partial sums + Richardson.
    Moderate precision; unreliable (flagged by callers) for stuples with (1,+1) blocks."""
    S = tuple(S)
    d = len(S)
    if d == 0: return mp.mpf(1)
    old = mp.mp.dps
    mp.mp.dps = dps + int(N*0.75) + 15
    M = 2*N + 8
    eps = [e for (s, e) in S]
    ss = [s for (s, e) in S]
    g = [mp.mpf(0)]*(M+1)
    acc = mp.mpf(0)
    ed, sd = eps[-1], ss[-1]
    for n in range(1, M+1):
        acc += (mp.mpf(ed)**n)/mp.mpf(n)**sd
        g[n] = acc
    for k in range(d-2, 0, -1):
        ek, sk = eps[k], ss[k]
        new = [mp.mpf(0)]*(M+1)
        acc = mp.mpf(0)
        for n in range(1, M+1):
            acc += (mp.mpf(ek)**n)*g[n-1]/mp.mpf(n)**sk
            new[n] = acc
        g = new
    e1, s1 = eps[0], ss[0]
    partials = []
    acc = mp.mpf(0)
    for n in range(1, M):
        if d == 1:
            term = (mp.mpf(e1)**n)/mp.mpf(n)**s1
        else:
            term = (mp.mpf(e1)**n)*g[n-1]/mp.mpf(n)**s1
        acc += term
        partials.append(acc)
    # Euler-average twice to suppress alternating sub-tails, then Richardson
    av = partials
    for _ in range(2):
        av = [(av[i]+av[i+1])/2 for i in range(len(av)-1)]
    val, maxc = mp.richardson(av)
    mp.mp.dps = old
    return +val

'''
src = src[:start] + new + src[end:]
p.write_text(src)
print("E3 v3")

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
    ((0,1,0,1), None),
    ((0,0,1,1), mp.pi**4/360),
    ((0,-1,-1), mp.zeta(3)/8),
    ((0,0,-1,-1), None),
    ((-1,-1,1), None),
    ((0,0,1), -mp.zeta(3)),
]
for w, ref in checks:
    S = rul1.word_to_stuple(w)
    v = (-1)**len(S)*rul1num.zeta_stuple_num(S, dps=45)
    line = "%s -> %s" % (w, mp.nstr(v,16))
    if ref is not None: line += "  err %s" % mp.nstr(abs(v-ref),3)
    print(line)
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=900)
print(r.stdout); print(r.stderr[-500:])

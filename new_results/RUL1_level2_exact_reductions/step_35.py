import pathlib
p = pathlib.Path('rul1num.py')
src = p.read_text()
start = src.index("def zeta_stuple_num")
end = src.index("def word_num")
new = '''def zeta_stuple_num(S, dps=40, N=160):
    """zeta(S) for convergent S: tail arrays + Richardson on the first 2N partial sums.
    Working precision raised to absorb Richardson cancellation."""
    S = tuple(S)
    d = len(S)
    if d == 0: return mp.mpf(1)
    old = mp.mp.dps
    mp.mp.dps = dps + int(N*0.75) + 15
    M = 2*N + 3
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
    val, maxc = mp.richardson(partials)
    mp.mp.dps = old
    return +val

'''
src = src[:start] + new + src[end:]
p.write_text(src)
print("rewritten")

import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import sys; sys.dont_write_bytecode=True
import mpmath as mp, time
import rul1num, rul1
mp.mp.dps = 50
checks = [
    ("z2",      ((2,1),), mp.pi**2/6),
    ("z3",      ((3,1),), mp.zeta(3)),
    ("z2bar",   ((2,-1),), -mp.pi**2/12),
    ("z21bar",  ((2,-1),(1,1)), mp.zeta(3)/8),
    ("z(1,1;-1,-1)", ((1,-1),(1,-1)), (mp.log(2)**2-mp.pi**2/6)/2),
    ("z(2,2)",  ((2,1),(2,1)), mp.pi**4/120),
    ("z(2,1)",  ((2,1),(1,1)), mp.zeta(3)),
    ("z(3,1)",  ((3,1),(1,1)), mp.pi**4/360),
    ("z(2,1,1)", ((2,1),(1,1),(1,1)), mp.pi**4/90),
]
for name, S, val in checks:
    t=time.time()
    got = rul1num.zeta_stuple_num(S, dps=45)
    print("E3 %-16s err %8s  (%.1fs)" % (name, mp.nstr(abs(got-val),3), time.time()-t))
rel = rul1.relation_DS((0,-1), (0,-1))
print("DS residual:", mp.nstr(abs(rul1num.combo_num(rel, dps=45)),3))
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=900)
print(r.stdout); print(r.stderr[-800:])

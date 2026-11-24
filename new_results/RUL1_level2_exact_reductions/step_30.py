import subprocess, pathlib
PY = '<project>/.venv/bin/python'
fix = r'''
import pathlib
p = pathlib.Path('rul1num.py')
src = p.read_text()
import re
# replace the whole zeta_stuple_num function
start = src.index("def zeta_stuple_num")
end = src.index("def word_num")
new = '''def zeta_stuple_num(S, dps=40, K=14, M=500):
    """zeta(S) for convergent S via tail arrays + Richardson on consecutive partial sums.
    Moderate precision (log-tail cases get fewer digits)."""
    old = mp.mp.dps
    mp.mp.dps = dps + 20
    S = tuple(S)
    d = len(S)
    if d == 0:
        mp.mp.dps = old; return mp.mpf(1)
    N = 2**K
    eps = [e for (s, e) in S]
    ss = [s for (s, e) in S]
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
    e1, s1 = eps[0], ss[0]
    from collections import deque
    win = deque(maxlen=M)
    acc = mp.mpf(0)
    for n in range(1, N+1):
        if d == 1:
            term = (mp.mpf(e1)**n)/mp.mpf(n)**s1
        else:
            term = (mp.mpf(e1)**n)*g[n-1]/mp.mpf(n)**s1
        acc += term
        win.append(acc)
    val, maxc = mp.richardson(list(win))
    mp.mp.dps = old
    return +val

'''
src = src[:start] + new + src[end:]
p.write_text(src)
print("rewritten zeta_stuple_num")
'''
print(subprocess.run([PY,'-c',fix], capture_output=True, text=True).stdout)
test = r'''
import sys; sys.dont_write_bytecode=True
import mpmath as mp
import rul1num, rul1
mp.mp.dps = 50
checks = [
    ('z2',      ((2,1),), mp.pi**2/6),
    ('z3',      ((3,1),), mp.zeta(3)),
    ('z2bar',   ((2,-1),), -mp.pi**2/12),
    ('z21bar',  ((2,-1),(1,1)), mp.zeta(3)/8),
    ('z(1,1;-1,-1)', ((1,-1),(1,-1)), (mp.log(2)**2-mp.pi**2/6)/2),
    ('z(2,2)',  ((2,1),(2,1)), mp.pi**4/120),
    ('z(2,1)',  ((2,1),(1,1)), mp.zeta(3)),
    ('z(3,1)',  ((3,1),(1,1)), mp.pi**4/360),
    ('z(2,1,1)', ((2,1),(1,1),(1,1)), mp.pi**4/90),
]
import time
for name, S, val in checks:
    t=time.time()
    got = rul1num.zeta_stuple_num(S, dps=45)
    print("E3 %-16s err %8s  (%.1fs)" % (name, mp.nstr(abs(got-val),3), time.time()-t))
rel = rul1.relation_DS((0,-1), (0,-1))
print("DS residual:", mp.nstr(abs(rul1num.combo_num(rel, dps=45)),3))
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=900)
print(r.stdout); print(r.stderr[-800:])

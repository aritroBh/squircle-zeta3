import pathlib
p = pathlib.Path('rul1num.py')
src = p.read_text()
start = src.index("def zeta_stuple_num")
end = src.index("def word_num")
new = '''def _neville_richardson(a):
    """Neville tableau for a_n = L + sum c_j/n^j (consecutive n = 1..M). Returns L."""
    R = list(a)
    M = len(R)
    for k in range(1, M):
        new = []
        for i in range(M-k):
            num = mp.mpf(i+1+k)**k * R[i+1] - mp.mpf(i+1)**k * R[i]
            den = mp.mpf(i+1+k)**k - mp.mpf(i+1)**k
            new.append(num/den)
        R = new
    return R[0]

def _euler_transform(b):
    """sum_{n>=1} (-1)^{n-1} b_n = sum_k D_k / 2^{k+1}, D_k = (Delta^k b)(1)."""
    D = list(b)
    tot = mp.mpf(0)
    pw = mp.mpf(1)/2
    k = 0
    while len(D) >= 1:
        tot += D[0]*pw
        pw /= 2
        D = [D[i]-D[i+1] for i in range(len(D)-1)]
        k += 1
    return tot

def zeta_stuple_num(S, dps=40, N=None):
    """zeta(S) for convergent S: tail arrays + own Neville/Euler acceleration."""
    S = tuple(S)
    d = len(S)
    if d == 0: return mp.mpf(1)
    old = mp.mp.dps
    if N is None: N = 2*int(dps) + 60
    mp.mp.dps = dps + 30
    e1, s1 = S[0]
    M = N + 5
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
    def b(n):
        if d == 1: return 1/mp.mpf(n)**s1
        return g[n-1]/mp.mpf(n)**s1
    if e1 == -1:
        val = -_euler_transform([b(n) for n in range(1, M)])
    else:
        val = _neville_richardson([mp.fsum([b(j) for j in range(1, n+1)]) for n in range(1, M)])
    mp.mp.dps = old
    return +val

'''
src = src[:start] + new + src[end:]
p.write_text(src)
print("E3 v4")

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
    ((-1,-1,1), -0.0947530042301788),
    ((0,0,1), -mp.zeta(3)),
    ((0,1), -mp.pi**2/6),
    ((0,1,1), mp.zeta(3)),
    ((0,1,1,1), mp.pi**4/90),
    ((-1,-1,-1), mp.log(2)**3/6),
    ((0,-1,0,1), -0.202935606320838410909250693101),
]
for w, ref in checks:
    S = rul1.word_to_stuple(w)
    v = (-1)**len(S)*rul1num.zeta_stuple_num(S, dps=45)
    line = "%s -> %s" % (w, mp.nstr(v,14))
    if ref is not None: line += "  err %s" % mp.nstr(abs(v-ref),3)
    print(line)
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=900)
print(r.stdout); print(r.stderr[-500:])

step = r'''
# step_03.py -- numeric validation of HPL evaluation and H <-> G sign conventions
import sys; sys.path.insert(0, '.')
import mpmath as mp
from hpl_num import hpl
mp.mp.dps = 60

def G_direct(letters, xx):
    """iterated integral G(letters; x) = int_0^x dt1/(t1-a1) int_0^{t1} ... via direct quadrature,
       only for convergent words (letters with a1 != 0 or inner vanishing)."""
    a = letters
    n = len(a)
    if n == 1:
        if a[0] == 0: return mp.log(xx)
        return mp.log((xx - a[0]) / (-a[0]))
    inner = lambda t: G_direct(a[1:], t)
    f = lambda t: inner(t) / (t - a[0])
    return mp.quad(f, [0, xx])

print("### HPL evaluator vs mpmath polylogs (Task 1(b) numeric)")
for X in (mp.mpf('0.3'), mp.mpf('0.866025403784')):
    e1 = hpl((0,0,0,1), X) - mp.polylog(4, X)
    e2 = hpl((0,0,1), X) - mp.polylog(3, X)
    e3 = hpl((0,1), X) - mp.polylog(2, X)
    e4 = hpl((0,-1), X) + mp.polylog(2, -X)
    print(f"  x={X}: Li4 {mp.nstr(e1,3)}  Li3 {mp.nstr(e2,3)}  Li2 {mp.nstr(e3,3)}  Li2(-x) {mp.nstr(e4,3)}")

print("### H_{1,0,1}(x) vs double sum Li_{1,2}(x) = sum_{k<l} x^l/(k l^2)")
for X in (mp.mpf('0.3'),):
    s = mp.nsum(lambda k, l: X**l / (k * l**2), [1, mp.inf], [lambda k: k + 1, mp.inf]) if False else None
    # direct truncated double sum
    acc = mp.mpf(0)
    N = 400
    for l in range(2, N):
        H1 = sum(mp.mpf(1)/k for k in range(1, l))
        acc += X**l * H1 / l**2
    tail = mp.mpf(0)  # x^400 tiny
    print("  H_{1,0,1}(0.3) =", mp.nstr(hpl((1,0,1), X), 20), "  sum =", mp.nstr(acc, 20),
          "  diff =", mp.nstr(hpl((1,0,1), X) - acc, 3))

print("### H <-> G sign convention: H_w = (-1)^{#1s} G_w (direct quadrature)")
for w in [(1,-1), (-1,1), (0,1,-1), (0,-1,1), (1,0,1), (1,0,-1), (0,0,-1), (0,1,1)]:
    X = mp.mpf('0.3')
    try:
        g = G_direct(w, X)
    except Exception as ex:
        print("  word", w, "quad failed:", ex); continue
    h = hpl(w, X)
    sgn = (-1)**sum(1 for t in w if t == 1)
    print(f"  word {w}: H={mp.nstr(h,15)}  (-1)^#1 G={mp.nstr(sgn*g,15)}  diff={mp.nstr(h-sgn*g,3)}")
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/step_03.py","w").write(step)
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_03.py"],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols")
print(r.stdout); print("STDERR:", (r.stderr or "")[-2000:])

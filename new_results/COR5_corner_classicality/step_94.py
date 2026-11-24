import sys, math
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
import mpmath as mp
from hpl_series import H_series

mp.mp.dps = 60
def f(a, t):
    return 1/t if a == 0 else (1/(1 - t) if a == 1 else 1/(1 + t))

def H_quad(w, x):
    if len(w) == 1:
        a = w[0]
        return mp.quad(lambda t: f(a, t), [0, x])
    a, rest = w[0], w[1:]
    return mp.quad(lambda t: f(a, t) * H_quad(rest, t), [0, x/4, x/2, 3*x/4, x])

c = mp.mpf(1)/2
for w in [(0,0,-1,1), (1,0,0,-1), (0,0,1,-1), (0,1,-1,1), (-1,1,-1,1)]:
    vq = H_quad(w, c)
    vs = H_series(w, c)
    print(f"H{w}: quad {mp.nstr(vq, 25)}")
    print(f"      series {mp.nstr(vs, 25)}  |diff| {mp.nstr(abs(vq-vs),3)}")

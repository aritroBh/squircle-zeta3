t = r'''
import sys; sys.path.insert(0, '.')
import mpmath as mp
from hpl_num import hpl
mp.mp.dps = 40
x = 1/mp.sqrt(2)

def G_rec(word, xx):
    word = tuple(word)
    if len(word) == 1:
        a = word[0]
        return mp.log(xx) if a == 0 else mp.log((xx - a)/(-a))
    a = word[0]
    return mp.quad(lambda t: G_rec(word[1:], t) / (t - a), [0, xx])

for w in [(1,-1,1,-1), (-1,1,-1,1)]:
    g = G_rec(w, x)
    h = hpl(w, x)
    sgn = (-1)**sum(1 for q in w if q == 1)
    print(w, " G(quad) =", mp.nstr(g, 25), " H(series) =", mp.nstr(h, 25),
          " sgn*G - H =", mp.nstr(sgn*g - h, 4))
# also the inner words
for w in [(-1,1,-1), (1,-1,1), (-1,1), (1,-1)]:
    g = G_rec(w, x)
    h = hpl(w, x)
    sgn = (-1)**sum(1 for q in w if q == 1)
    print(w, " G(quad) =", mp.nstr(g, 25), " H(series) =", mp.nstr(h, 25),
          " sgn*G - H =", mp.nstr(sgn*g - h, 4))
'''
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=900)
print(r.stdout); print("ERR:", (r.stderr or "")[-500:])

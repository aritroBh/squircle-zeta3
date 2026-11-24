t = r'''
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
import mpmath as mp, json, re
from hpl_num import hpl
mp.mp.dps = 40
data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
W4 = sorted(set(parse_term(s)[1] for s, c in data if len(parse_term(s)[1]) == 4))
def G_rec(word, xx):
    word = tuple(word)
    if len(word) == 1:
        a = word[0]
        return mp.log(xx) if a == 0 else mp.log((xx - a)/(-a))
    a = word[0]
    return mp.quad(lambda t: G_rec(word[1:], t) / (t - a), [0, xx])
x = mp.mpf(1)/2
bad = 0
for w in W4:
    g = G_rec(w, x)
    h = hpl(w, x)
    sgn = (-1)**sum(1 for q in w if q == 1)
    err = abs(h - sgn*g)
    if err > mp.mpf('1e-25'):
        bad += 1
        print("SIGN MISMATCH", w, mp.nstr(h, 12), mp.nstr(sgn*g, 12), mp.nstr(err, 3))
print("words with sign mismatch:", bad, "of", len(W4))
'''
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=1200)
print(r.stdout); print("ERR:", (r.stderr or "")[-500:])
